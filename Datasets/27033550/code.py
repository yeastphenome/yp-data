#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 27033550
paper_name = 'yang_nystrom_2016' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[15]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, names=['orfs','values'], sep='\t')


# In[16]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[17]:


original_data['orfs'] = original_data['orfs'].astype(str)


# In[18]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[19]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[20]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[21]:


original_data['values'] = original_data['values'].astype(float)


# In[22]:


# Log-transform fold enrichments so that the mode (all non-hit strains) will be 0, by default
original_data['data'] = np.log2(original_data['values'])


# In[23]:


original_data.set_index('orfs', inplace=True)
original_data.index.name='orf'


# In[24]:


original_data = original_data[['data']].copy()


# In[25]:


original_data = original_data.groupby(original_data.index).mean()


# In[26]:


original_data.head()


# # Prepare the final dataset

# In[27]:


dataset_ids = [16413]


# In[28]:


datasets = datasets.reindex(index=dataset_ids)


# In[29]:


data = original_data.copy()


# In[30]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[32]:


data.head()


# ## Subset to the genes currently in SGD

# In[33]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[34]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[35]:


data.head()


# In[36]:


data.shape


# # Normalize

# In[37]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[38]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[39]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[40]:


data_all.head()


# # Print out

# In[41]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[42]:


from IO.save_data_to_db3 import *


# In[43]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




