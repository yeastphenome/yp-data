#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 14661025
paper_name = 'parsons_boone_2004' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[11]:


original_data = pd.read_excel('raw_data/nbt919-S2.xlsx', sheet_name='Sheet1', skiprows=6)


# In[12]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[13]:


original_data.head()


# In[14]:


original_data['orf'] = original_data['ORF'].astype(str)


# In[15]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[16]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[17]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[18]:


original_data.set_index('orf', inplace=True)


# In[19]:


original_data.drop(columns=['ORF','GENE'], inplace=True)


# In[20]:


original_data = original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[21]:


original_data[original_data.isnull()] = 0


# In[22]:


original_data = original_data.groupby(original_data.index).mean()


# In[23]:


original_data.shape


# In[24]:


# Flip the sign of the values, such that negative = sensitive, positive = resistant
original_data = -original_data


# # Load datasets

# In[27]:


dt = pd.read_csv('extras/datasets.txt', sep='\t', header=None)


# In[28]:


dt.head()


# In[29]:


dt.set_index(1, inplace=True)


# In[30]:


dt = dt.reindex(index=original_data.columns.values)


# In[32]:


dataset_ids = dt[0].values


# # Prepare the final dataset

# In[33]:


data = original_data.copy()


# In[34]:


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


data_norm = normalize_phenotypic_scores(data, has_tested=False)


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




