#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 23689276
paper_name = 'bowie_fyles_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[20]:


original_data1 = pd.read_excel('raw_data/c3ob40593a.xlsx', sheet_name='Table S1', skiprows=2)
original_data2 = pd.read_excel('raw_data/c3ob40593a.xlsx', sheet_name='Table S2', skiprows=2)


# In[21]:


print('Original data dimensions: %d x %d' % (original_data1.shape))


# In[22]:


original_data1['Orf'] = original_data1['Orf'].astype(str)
original_data2['Orf'] = original_data2['Orf'].astype(str)


# In[23]:


# Eliminate all white spaces & capitalize
original_data1['Orf'] = clean_orf(original_data1['Orf'])
original_data2['Orf'] = clean_orf(original_data2['Orf'])


# In[24]:


# Translate to ORFs 
original_data1['Orf'] = translate_sc(original_data1['Orf'], to='orf')
original_data2['Orf'] = translate_sc(original_data2['Orf'], to='orf')


# In[25]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['Orf'])
print(original_data1.loc[~t,])


# In[26]:


original_data1 = original_data1.loc[t,:]


# In[27]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['Orf'])
print(original_data2.loc[~t,])


# In[28]:


original_data2 = original_data2.loc[t,:]


# In[29]:


original_data1.set_index('Orf', inplace=True)
original_data2.set_index('Orf', inplace=True)


# In[30]:


original_data1.index.name='orf'
original_data2.index.name='orf'


# In[32]:


original_data1['data'] = pd.to_numeric(original_data1['Ratio'], errors='coerce')
original_data2['data'] = pd.to_numeric(original_data2['Ratio'], errors='coerce')


# In[33]:


original_data1 = original_data1[['data']].copy()
original_data2 = original_data2[['data']].copy()


# In[34]:


original_data1 = original_data1.groupby(original_data1.index).mean()
original_data2 = original_data2.groupby(original_data2.index).mean()


# In[35]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[36]:


original_data.head()


# In[37]:


original_data.shape


# # Prepare the final dataset

# In[38]:


data = original_data.copy()


# In[39]:


dataset_ids = [16557, 16558]


# In[40]:


datasets = datasets.reindex(index=dataset_ids)


# In[41]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[42]:


data.head()


# ## Subset to the genes currently in SGD

# In[43]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[44]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[45]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[46]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[47]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[48]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[49]:


from IO.save_data_to_db3 import *


# In[50]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




