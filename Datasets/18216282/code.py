#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 18216282
paper_name = 'schluter_conibear_2008' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[38]:


original_data = pd.read_excel('raw_data/DataSet1.xlsx', sheet_name='Table S1', skiprows=2)


# In[39]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[40]:


original_data.head()


# In[41]:


original_data['orf'] = original_data['ORF'].astype(str)


# In[42]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[43]:


original_data.loc[original_data['orf']=='YLR287-A','orf'] = 'YLR287C-A'


# In[44]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[45]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[46]:


original_data.set_index('orf', inplace=True)


# In[47]:


original_data = original_data[['MATa 1','MATa 2','MATa 1.1','MATa 2.1','diploid 1','diploid 2']].copy()


# In[48]:


original_data = -original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[49]:


original_data = original_data.groupby(original_data.index).mean()


# In[50]:


original_data.shape


# In[51]:


original_data.columns = [16030, 16030, 16031, 16031, 16032, 16032]


# In[52]:


original_data = original_data.T
original_data = original_data.groupby(original_data.index).mean()
original_data = original_data.T


# In[53]:


original_data.shape


# # Prepare the final dataset

# In[54]:


data = original_data.copy()


# In[55]:


dataset_ids = original_data.columns.values
datasets = datasets.reindex(index=dataset_ids)


# In[56]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[57]:


data.head()


# ## Subset to the genes currently in SGD

# In[58]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[59]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[60]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[61]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[62]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[63]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[64]:


from IO.save_data_to_db3 import *


# In[65]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




