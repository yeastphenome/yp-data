#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[40]:


paper_pmid = 30647105
paper_name = 'alhoch_tang_2019' 


# In[41]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[42]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[43]:


original_data = pd.read_excel('raw_data/S6 SC genomic screen BHA and BPA TN1.xlsx', sheet_name='Sc BPA and BHA genomic', skiprows=1)


# In[44]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[45]:


original_data['ORF name'] = original_data['ORF name'].astype(str)


# In[46]:


# Eliminate all white spaces & capitalize
original_data['ORF name'] = clean_orf(original_data['ORF name'])


# In[47]:


# Translate to ORFs 
original_data['ORF name'] = translate_sc(original_data['ORF name'], to='orf')


# In[48]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF name'])
print(original_data.loc[~t,])


# In[49]:


original_data = original_data.loc[t,]


# In[50]:


original_data.set_index('ORF name', inplace=True)
original_data.index.name='orf'


# In[51]:


original_data[original_data.notnull()] = -1


# In[52]:


original_data[original_data.isnull()] = 0


# In[53]:


original_data.sum(axis=0)


# In[54]:


original_data = original_data.astype(float)


# In[55]:


original_data = original_data.groupby(original_data.index).mean()


# # Prepare the final dataset

# In[56]:


data = original_data.copy()


# In[57]:


dataset_ids = [16600, 16599]


# In[58]:


datasets = datasets.reindex(index=dataset_ids)


# In[59]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[60]:


data.head()


# ## Subset to the genes currently in SGD

# In[61]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[63]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[64]:


data.head()


# In[65]:


data.shape


# # Normalize

# In[66]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[67]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[68]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[69]:


data_all.head()


# # Print out

# In[70]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[71]:


from IO.save_data_to_db3 import *


# In[72]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




