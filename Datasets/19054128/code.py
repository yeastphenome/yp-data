#!/usr/bin/env python
# coding: utf-8

# In[5]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[6]:


paper_pmid = 19054128
paper_name = 'yoshikawa_shimizu_2009' 


# In[7]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[8]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[12]:


original_data = pd.read_excel('raw_data/FYR_456_sm_tableS1.xlsx', sheet_name='data', skiprows=1)


# In[13]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[14]:


original_data.head()


# In[15]:


original_data['orf'] = original_data['Name'].astype(str)


# In[16]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[17]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[18]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[19]:


original_data.set_index('orf', inplace=True)


# In[21]:


data_cols = [c for c in original_data.columns if 'specific growth rate' in c]


# In[22]:


original_data = original_data[data_cols].copy()


# In[23]:


original_data = original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[24]:


original_data = original_data.groupby(original_data.index).mean()


# In[25]:


original_data.shape


# In[27]:


dataset_ids = [100, 100, 523, 523, 4, 4, 5, 5]
original_data.columns = dataset_ids


# In[28]:


original_data = original_data.T
original_data = original_data.groupby(original_data.index).mean()
original_data = original_data.T


# In[29]:


original_data.shape


# In[31]:


# Normalize by the untreated sample
original_data = original_data.div(original_data[100], axis=0)
original_data.head()


# In[32]:


original_data.drop(columns=[100], inplace=True)


# # Prepare the final dataset

# In[33]:


data = original_data.copy()


# In[34]:


dataset_ids = original_data.columns.values
datasets = datasets.reindex(index=dataset_ids)


# In[35]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[36]:


data.head()


# ## Subset to the genes currently in SGD

# In[37]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[38]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[39]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[40]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[41]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[42]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[43]:


from IO.save_data_to_db3 import *


# In[44]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




