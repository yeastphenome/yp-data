#!/usr/bin/env python
# coding: utf-8

# In[12]:


get_ipython().run_line_magic('run', '../yp_utils.py')
import re


# # Initial setup

# In[2]:


paper_pmid = 31451498
paper_name = 'parisi_bleackley_2019' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/AAC.01097-19-sd002.xlsx', sheet_name='Log2 ratios')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[40]:


original_data['orf'] = original_data['BARCODE'].astype(str)


# In[41]:


strs = ['-1t','-2t','-1d','-2d','-1l','-1','-2','-3','-t','-d']


# In[42]:


# Remove "-1", "-2" and "-1t" from orfs
for s in strs:
    original_data['orf'] = original_data['orf'].apply(lambda x: x.replace(s,''))


# In[43]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[44]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'].values, to='orf')


# In[45]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[46]:


original_data.set_index('orf', inplace=True)


# In[47]:


original_data = original_data[['DmAMP1', 'NaD1', 'NbD6', 'SBI6']].copy()
original_data = original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[48]:


original_data = original_data.groupby(original_data.index).mean()


# In[49]:


original_data.shape


# In[61]:


num_vals = (np.abs(original_data)>0).sum(axis=1)


# In[64]:


original_data = original_data.loc[num_vals > 0,:]


# In[65]:


original_data.shape


# # Prepare the final dataset

# In[66]:


data = original_data.copy()


# In[67]:


dataset_ids = [21834,21833,21835,21836]
datasets = datasets.reindex(index=dataset_ids)


# In[68]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[69]:


data.head()


# ## Subset to the genes currently in SGD

# In[70]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[71]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[72]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[73]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[74]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[75]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[76]:


from IO.save_data_to_db3 import *


# In[77]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




