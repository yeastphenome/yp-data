#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 111111111
paper_name = 'firstauthor_lastauthor_year' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, names=['genes'])


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[ ]:


original_data.head()


# In[7]:


original_data['orf'] = original_data['orf'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[9]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[12]:


original_data['data'] = -1


# In[13]:


original_data.set_index('orf', inplace=True)


# In[ ]:


original_data = original_data[['data']].copy()


# In[ ]:


original_data = original_data.groupby(original_data.index).mean()


# In[ ]:


original_data.shape


# # Load & process tested strains

# In[14]:


tested = pd.read_excel('raw_data/Mat_a_obs_v4 0.xlsx', sheet_name='DATA')


# In[15]:


tested.head()


# In[16]:


tested['orf'] = tested['ORF name'].astype(str)


# In[ ]:


tested['orf'] = clean_orf(tested['orf'])


# In[ ]:


tested.loc[tested['orf']=='YLR287-A','orf'] = 'YLR287C-A'


# In[ ]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[ ]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[ ]:


tested = tested.loc[t,:]


# In[ ]:


tested_orfs = tested['orf'].unique()


# In[ ]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[ ]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[ ]:


data = original_data.copy()


# In[17]:


dataset_ids = [16461]
datasets = datasets.reindex(index=dataset_ids)


# In[ ]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[26]:


data.head()


# ## Subset to the genes currently in SGD

# In[ ]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[ ]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[ ]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[ ]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[ ]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[27]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[2]:


from IO.save_data_to_db3 import *


# In[29]:


save_data_to_db(data_all, paper_pmid)

