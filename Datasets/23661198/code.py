#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 23661198
paper_name = 'vahey_voldman_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data1 = pd.read_excel('raw_data/c3lc50162k.xlsx', sheet_name='Sheet1', skiprows=1)
original_data2 = pd.read_excel('raw_data/c3lc50162k_2.xlsx', sheet_name='Sheet1', skiprows=1)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[7]:


original_data2.head()


# In[8]:


original_data1['orf'] = original_data1['Unnamed: 0'].astype(str)
original_data2['orf'] = original_data2['Unnamed: 0'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data1['orf'] = clean_orf(original_data1['orf'])
original_data2['orf'] = clean_orf(original_data2['orf'])


# In[12]:


# Translate to ORFs 
original_data1['orf'] = translate_sc(original_data1['orf'], to='orf')
original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')


# In[13]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orf'])
print(original_data1.loc[~t,])


# In[14]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[15]:


original_data1['data'] = original_data1['y (avg.)'].astype(float)
original_data2['data'] = original_data2['y (avg.)'].astype(float)


# In[16]:


original_data1.set_index('orf', inplace=True)
original_data2.set_index('orf', inplace=True)


# In[17]:


original_data1 = original_data1[['data']].copy()
original_data2 = original_data2[['data']].copy()


# In[18]:


original_data1 = original_data1.groupby(original_data1.index).mean()
original_data2 = original_data2.groupby(original_data2.index).mean()


# In[20]:


original_data = original_data1.join(original_data2, lsuffix='_1', rsuffix='_2', how='outer')


# In[21]:


original_data.head()


# # Prepare the final dataset

# In[22]:


data = original_data.copy()


# In[23]:


dataset_ids = [95, 16597]


# In[24]:


datasets = datasets.reindex(index=dataset_ids)


# In[25]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[26]:


data.head()


# ## Subset to the genes currently in SGD

# In[27]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[28]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# In[29]:


data.shape


# # Normalize

# In[30]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[31]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[32]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[33]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[34]:


from IO.save_data_to_db3 import *


# In[35]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




