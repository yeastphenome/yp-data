#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 24142923
paper_name = 'jarolim_dawes_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data1 = pd.read_excel('raw_data/FileS1.xlsx', sheet_name='Exponential_Resistant', skiprows=2)
original_data2 = pd.read_excel('raw_data/FileS1.xlsx', sheet_name='Exponential_Sensitive', skiprows=2)
original_data3 = pd.read_excel('raw_data/FileS1.xlsx', sheet_name='Stationary_Resistant', skiprows=2)
original_data4 = pd.read_excel('raw_data/FileS1.xlsx', sheet_name='Stationary_Sensitive', skiprows=2)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))
print('Original data dimensions: %d x %d' % (original_data3.shape))
print('Original data dimensions: %d x %d' % (original_data4.shape))


# In[7]:


orf_col = 'Systematic Name'
original_data1['orfs'] = original_data1['SystematicName'].astype(str)
original_data2['orfs'] = original_data2[orf_col].astype(str)
original_data3['orfs'] = original_data3[orf_col].astype(str)
original_data4['orfs'] = original_data4[orf_col].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data1['orfs'] = clean_orf(original_data1['orfs'])
original_data2['orfs'] = clean_orf(original_data2['orfs'])
original_data3['orfs'] = clean_orf(original_data3['orfs'])
original_data4['orfs'] = clean_orf(original_data4['orfs'])


# In[9]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['orfs'], to='orf')
original_data2['orfs'] = translate_sc(original_data2['orfs'], to='orf')
original_data3['orfs'] = translate_sc(original_data3['orfs'], to='orf')
original_data4['orfs'] = translate_sc(original_data4['orfs'], to='orf')


# In[10]:


# Fix typos
original_data3.loc[original_data3['orfs']=='YML095-A','orfs'] = 'YML095C-A'


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[13]:


# Make sure everything translated ok
t = looks_like_orf(original_data3['orfs'])
print(original_data3.loc[~t,])


# In[14]:


# Make sure everything translated ok
t = looks_like_orf(original_data4['orfs'])
print(original_data4.loc[~t,])


# In[15]:


original_data1.head()


# In[16]:


original_data1['data'] = original_data1['Rating']
original_data2['data'] = original_data2['Rating']
original_data3['data'] = original_data3['Rating']
original_data4['data'] = original_data4['Rating']


# In[17]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)
original_data3.set_index('orfs', inplace=True)
original_data4.set_index('orfs', inplace=True)


# In[18]:


original_data1.index.name = 'orf'
original_data2.index.name = 'orf'
original_data3.index.name = 'orf'
original_data4.index.name = 'orf'


# In[20]:


original_data1 = original_data1[['data']]
original_data2 = original_data2[['data']]
original_data3 = original_data3[['data']]
original_data4 = original_data4[['data']]


# In[21]:


original_data1 = original_data1.groupby(original_data1.index).mean()
original_data2 = original_data2.groupby(original_data2.index).mean()
original_data3 = original_data3.groupby(original_data3.index).mean()
original_data4 = original_data4.groupby(original_data4.index).mean()


# In[22]:


data1 = original_data1[['data']].join(original_data2[['data']], how='outer', lsuffix='_r', rsuffix='_s')


# In[23]:


data1['data'] = data1[['data_s','data_r']].mean(axis=1)


# In[24]:


data2 = original_data3[['data']].join(original_data4[['data']], how='outer', lsuffix='_r', rsuffix='_s')


# In[25]:


data2['data'] = data2[['data_s','data_r']].mean(axis=1)


# In[26]:


data = data1[['data']].join(data2[['data']], how='outer', lsuffix='_exp', rsuffix='_stn')


# In[27]:


data[data.isnull()] = 0


# In[28]:


data.head()


# # Prepare the final dataset

# In[29]:


dataset_ids = [16624,16538]
datasets = datasets.reindex(index=dataset_ids)


# In[30]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[31]:


data.head()


# ## Subset to the genes currently in SGD

# In[32]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[33]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# In[34]:


data.shape


# # Normalize

# In[35]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[36]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[37]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[38]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[39]:


from IO.save_data_to_db3 import *


# In[40]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




