#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[56]:


paper_pmid = 23959800
paper_name = 'papic_rapaport_2013' 


# In[57]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[58]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[44]:


original_data = pd.read_csv('raw_data/hits.txt', sep='\t')


# In[45]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[46]:


original_data.head()


# In[47]:


original_data['genes'] = original_data['Gene'].astype(str)


# In[48]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[49]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['genes'], to='orf')


# In[50]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[51]:


original_data.set_index('orf', inplace=True)


# In[52]:


original_data['data'] = -original_data[['Value 1','Value 2']].mean(axis=1)


# In[53]:


original_data = original_data[['data']].copy()


# In[54]:


original_data = original_data.groupby(original_data.index).mean()


# In[55]:


original_data.shape


# # Prepare the final dataset

# In[59]:


data = original_data.copy()


# In[60]:


dataset_ids = [11867]
datasets = datasets.reindex(index=dataset_ids)


# In[61]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[62]:


data.head()


# ## Subset to the genes currently in SGD

# In[63]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[64]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[65]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[66]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[67]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[68]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[69]:


from IO.save_data_to_db3 import *


# In[70]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




