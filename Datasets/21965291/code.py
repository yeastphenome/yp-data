#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 21965291
paper_name = 'zakrzewska_smits_2011' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[28]:


original_data = pd.read_excel('raw_data/mc-E10-08-0721-s06.xlsx', sheet_name='growth rates')


# In[29]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[30]:


original_data.head()


# In[31]:


original_data['orf'] = original_data['rc>0'].astype(str)


# In[32]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[33]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[34]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[35]:


original_data.set_index('orf', inplace=True)


# In[36]:


original_data = original_data[['mu 30','mu 38']].copy()


# In[37]:


for c in original_data.columns:
    original_data[c] = pd.to_numeric(original_data[c], errors='coerce')


# In[38]:


original_data = original_data.groupby(original_data.index).mean()


# In[39]:


original_data.shape


# # Dataset 2

# In[40]:


original_data2 = pd.read_excel('raw_data/mc-E10-08-0721-s06.xlsx', sheet_name='survival % 95% CI')


# In[41]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[42]:


original_data2.head()


# In[43]:


original_data2['orf'] = original_data2['rc>0'].astype(str)


# In[44]:


# Eliminate all white spaces & capitalize
original_data2['orf'] = clean_orf(original_data2['orf'])


# In[45]:


# Translate to ORFs 
original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')


# In[46]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[47]:


original_data2.set_index('orf', inplace=True)


# In[48]:


original_data2 = original_data2[['30oC oxi','30oC acid','30oC heat','38oC oxi','38oC acid','38oC heat']].copy()


# In[50]:


for c in original_data2.columns:
    original_data2[c] = pd.to_numeric(original_data2[c], errors='coerce')


# In[51]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[52]:


original_data2.shape


# # Merge

# In[53]:


original_data = original_data.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[54]:


original_data.head()


# # Prepare the final dataset

# In[55]:


data = original_data.copy()


# In[56]:


dataset_ids = list(np.arange(16128,16136))
datasets = datasets.reindex(index=dataset_ids)


# In[57]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[58]:


data.head()


# ## Subset to the genes currently in SGD

# In[59]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[60]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[61]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[62]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[63]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[64]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[65]:


from IO.save_data_to_db3 import *


# In[66]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




