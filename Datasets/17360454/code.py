#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 17360454
paper_name = 'yuen_spencer_2007' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


sheets = ['SI Table 3','SI Table 4']


# In[7]:


data_cols = ['CTF','BiM','ALF']


# In[9]:


original_data_list = []
for s in sheets:
    original_data = pd.read_excel('raw_data/10642Tables_3-9.xlsx', sheet_name=s)
    print('Original data dimensions: %d x %d' % (original_data.shape))
#     print(original_data.head())
    original_data['orf'] = original_data['ORF name'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data.set_index('orf', inplace=True)
    original_data = original_data[data_cols].copy()
    original_data = original_data.apply(pd.to_numeric, axis=1, errors='coerce')
    original_data = original_data.groupby(original_data.index).mean()
    original_data = original_data - 1
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[10]:


original_data = pd.concat(original_data_list, axis=0)


# In[12]:


original_data.shape


# In[13]:


original_data.head()


# # Prepare the final dataset

# In[14]:


data = original_data.copy()


# In[15]:


dataset_ids = [4931, 4932, 4933]
datasets = datasets.reindex(index=dataset_ids)


# In[16]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[17]:


data.head()


# ## Subset to the genes currently in SGD

# In[18]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[19]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[20]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[21]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[22]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[23]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[24]:


from IO.save_data_to_db3 import *


# In[25]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:



