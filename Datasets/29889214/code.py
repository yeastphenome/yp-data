#!/usr/bin/env python
# coding: utf-8

# In[11]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[12]:


paper_pmid = 29889214
paper_name = 'wong_khalil_2018' 


# In[13]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[14]:


datasets.set_index('dataset_id', inplace=True)


# In[25]:


datasets.shape


# # Load & process the data

# In[15]:


original_data = pd.read_csv('raw_data/YKO_HS_fitness.csv', sep=',')


# In[16]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[17]:


original_data.head()


# In[20]:


original_data.drop(columns=['Temperature Magnitude'], inplace=True)


# In[22]:


original_data.drop(index=[0,1,2,3], inplace=True)


# In[26]:


original_data['orf'] = original_data['Unnamed: 0'].astype(str)


# In[27]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[28]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[29]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[30]:


original_data.set_index('orf', inplace=True)


# In[31]:


original_data.drop(columns=['Unnamed: 0'], inplace=True)


# In[32]:


original_data.columns = datasets.index.values


# In[33]:


original_data = original_data.apply(pd.to_numeric, axis=1)


# In[34]:


original_data = original_data.groupby(original_data.index).mean()


# In[35]:


original_data.shape


# # Prepare the final dataset

# In[36]:


data = original_data.copy()


# In[37]:


dataset_ids = original_data.columns.values
datasets = datasets.reindex(index=dataset_ids)


# In[38]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[39]:


data.head()


# ## Subset to the genes currently in SGD

# In[40]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[41]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[42]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[43]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[44]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[45]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[46]:


from IO.save_data_to_db3 import *


# In[47]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




