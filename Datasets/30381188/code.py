#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[3]:


paper_pmid = 30381188
paper_name = 'ruta_farcasanu_2018' 


# In[4]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[5]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[11]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, names=['orf'])


# In[12]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[13]:


original_data['orf'] = original_data['orf'].astype(str)


# In[14]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[15]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[17]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[18]:


original_data['data'] = -1


# In[19]:


original_data.set_index('orf', inplace=True)


# In[20]:


original_data = original_data.groupby(original_data.index).mean()


# In[21]:


original_data.shape


# # Prepare the final dataset

# In[22]:


data = original_data[['data']].copy()


# In[23]:


dataset_ids = [16209]
datasets = datasets.reindex(index=dataset_ids)


# In[24]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[25]:


data.head()


# ## Subset to the genes currently in SGD

# In[26]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[27]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[28]:


data.head()


# # Normalize

# In[30]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[31]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[32]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[33]:


data_all.head()


# # Print out

# In[34]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[35]:


from IO.save_data_to_db3 import *


# In[36]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




