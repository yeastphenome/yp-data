#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 19619494
paper_name = 'okamoto_ohsumi_2009' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# ### Part1

# In[5]:


original_data1 = pd.read_csv('raw_data/Table1.txt', header=None, names=['genes','data'], sep='\t')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))


# In[7]:


original_data1['genes'] = original_data1['genes'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data1['genes'] = clean_genename(original_data1['genes'])


# In[9]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['genes'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[11]:


data_dict = {'+++': 0, '++': -1, '+': -2, '-': -3}


# In[12]:


original_data1['data_score'] = original_data1['data'].apply(lambda x: data_dict[x])


# In[13]:


original_data1['data_score'] = original_data1['data_score'].astype(int)


# In[14]:


original_data1.set_index('orfs', inplace=True)


# ### Part2

# In[15]:


original_data2 = pd.read_csv('raw_data/TableS1.txt', header=None, names=['genes','data'], sep='\t')


# In[16]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[17]:


original_data2['genes'] = original_data2['genes'].astype(str)


# In[18]:


# Eliminate all white spaces & capitalize
original_data2['genes'] = clean_genename(original_data2['genes'])


# In[19]:


# Translate to ORFs 
original_data2['orfs'] = translate_sc(original_data2['genes'], to='orf')


# In[20]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[21]:


data_dict = {'+++': 0, '++': -1, '+': -2, '-': -3}


# In[22]:


original_data2['data_score'] = original_data2['data'].apply(lambda x: data_dict[x])


# In[23]:


original_data2['data_score'] = original_data2['data_score'].astype(int)


# In[24]:


original_data2.set_index('orfs', inplace=True)


# In[25]:


original_data2.head()


# # Merge

# In[26]:


original_data1.index.name = 'orf'
original_data2.index.name = 'orf'

original_data1 = original_data1[['data_score']].copy()
original_data2 = original_data2[['data_score']].copy()


# In[27]:


original_data1 = original_data1.groupby(original_data1.index).mean()
original_data2 = original_data2.groupby(original_data2.index).mean()


# In[28]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_s1')


# In[30]:


original_data['data'] = original_data.mean(axis=1)


# In[31]:


original_data.head()


# In[32]:


original_data = original_data[['data']].copy()


# # Prepare the final dataset

# In[33]:


data = original_data.copy()


# In[34]:


dataset_ids = [16606]
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

# In[40]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[41]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[42]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[43]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[44]:


from IO.save_data_to_db3 import *


# In[45]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




