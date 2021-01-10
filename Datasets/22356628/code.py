#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 22356628
paper_name = 'davey_griffith_2012' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data1 = pd.read_csv('raw_data/emi_2705_sm_tables1.txt', header=None, names=['orfs'])
original_data2 = pd.read_csv('raw_data/emi_2705_sm_tables2.txt', header=None, names=['orfs'])


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[7]:


original_data1['orfs'] = original_data1['orfs'].astype(str)
original_data2['orfs'] = original_data2['orfs'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data1['orfs'] = clean_orf(original_data1['orfs'])
original_data2['orfs'] = clean_orf(original_data2['orfs'])


# In[9]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['orfs'], to='orf')
original_data2['orfs'] = translate_sc(original_data2['orfs'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[12]:


original_data1['data'] = 1
original_data2['data'] = -1


# In[13]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)
original_data1.index.name='orf'
original_data2.index.name='orf'


# In[16]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[17]:


original_data[original_data.isnull()] = 0


# # Prepare the final dataset

# In[21]:


data = original_data.copy()


# In[22]:


dataset_ids = [16517, 16513]
datasets = datasets.reindex(index=dataset_ids)


# In[23]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[24]:


data.head()


# ## Subset to the genes currently in SGD

# In[25]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[26]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[28]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[29]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[30]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[31]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[32]:


from IO.save_data_to_db3 import *


# In[33]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




