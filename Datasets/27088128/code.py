#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 27088128
paper_name = 'koselny_krysan_2016' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[23]:


original_data1 = pd.read_excel("raw_data/deliv_AltComparisons - Lacey's data - Sorted by Pvalue.xlsx", sheet_name='DowntagFC', skiprows=1)
original_data2 = pd.read_excel("raw_data/deliv_AltComparisons - Lacey's data - Sorted by Pvalue.xlsx", sheet_name='UptagFC', skiprows=1)


# In[24]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[25]:


original_data1.head()


# In[26]:


original_data1['orf'] = original_data1['orf_name'].astype(str)
original_data2['orf'] = original_data2['orf_name'].astype(str)


# In[27]:


# Eliminate all white spaces & capitalize
original_data1['orf'] = clean_orf(original_data1['orf'])
original_data2['orf'] = clean_orf(original_data2['orf'])


# In[28]:


# Translate to ORFs 
original_data1['orf'] = translate_sc(original_data1['orf'], to='orf')
original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')


# In[29]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orf'])
print(original_data1.loc[~t,])


# In[30]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[31]:


original_data1.loc[original_data1['orf'] == 'YBR160W_AS','orf'] = 'YBR160W'
original_data2.loc[original_data2['orf'] == 'YBR160W_AS','orf'] = 'YBR160W'


# In[32]:


for c in ['A','B','C','D']:
    original_data1[c] = pd.to_numeric(original_data1[c], errors='coerce')
    original_data2[c] = pd.to_numeric(original_data2[c], errors='coerce')


# In[33]:


original_data1 = original_data1.groupby(original_data1['orf']).mean()
original_data2 = original_data2.groupby(original_data2['orf']).mean()


# In[34]:


print(original_data1.shape)
print(original_data2.shape)


# In[37]:


original_data = original_data1[['A','B','C','D']].join(original_data2[['A','B','C','D']], lsuffix='_down', rsuffix='_up')


# In[38]:


original_data.head()


# In[39]:


original_data['data'] = original_data.mean(axis=1)


# In[41]:


original_data.head()


# # Prepare the final dataset

# In[54]:


data = original_data[['data']].copy()


# In[55]:


dataset_ids = [22077]
datasets = datasets.reindex(index=dataset_ids)


# In[56]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[57]:


data.head()


# ## Subset to the genes currently in SGD

# In[58]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[59]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[60]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[61]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[62]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[63]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[64]:


from IO.save_data_to_db3 import *


# In[65]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:



