#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 31511699
paper_name = 'puddu_jackson_2019' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[7]:


original_data = pd.read_excel('raw_data/41586_2019_1549_MOESM3_ESM.xlsx', 
                              sheet_name='SupplementaryTable3', skiprows=33)


# In[8]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[9]:


original_data.head()


# In[10]:


hit_strains_ids = original_data.iloc[:,12]


# In[16]:


# Extract gene names
gene_names = [x.split('_')[1] for x in hit_strains_ids.values]


# In[18]:


original_data['genes'] = gene_names
original_data['genes'] = original_data['genes'].astype(str)


# In[19]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[20]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[24]:


# Fix a few problems manually
rename_map = {'FLO8': 'YER109C','AAD6':'YFL056C','SDL1':'YIL167W','HXT12':'YIL170W','SDC25':'YLL016W','CRS5':'YOR031W'}
original_data['orfs'] = original_data['orfs'].apply(lambda x: rename_map[x] if x in rename_map.keys() else x)


# In[25]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[26]:


original_data = original_data.loc[t,:]


# In[27]:


original_data.set_index('orfs', inplace=True)
original_data.index.name='orf'


# In[39]:


columns = np.concatenate((original_data.columns.values[1:11], original_data.columns.values[13:29], original_data.columns.values[[30]]))


# In[42]:


original_data = original_data.loc[:,columns]


# In[43]:


original_data = original_data.groupby(original_data.index).mean()


# In[44]:


original_data.shape


# # Prepare the final dataset

# In[46]:


data = original_data.copy()


# In[48]:


dataset_ids = np.arange(16322,16349)
datasets = datasets.reindex(index=dataset_ids)


# In[51]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[52]:


data.head()


# ## Subset to the genes currently in SGD

# In[53]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[54]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[55]:


data.head()


# # Normalize

# In[56]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[57]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[58]:


data_norm[data.isnull()] = np.nan


# In[59]:


data_all = data.join(data_norm)


# In[60]:


data_all.head()


# # Print out

# In[61]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[62]:


from IO.save_data_to_db3 import *


# In[63]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




