#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 20657825
paper_name = 'fabrizio_longo_2010' 


# In[46]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[47]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[27]:


original_data = pd.read_excel('raw_data/journal.pgen.1001024.s006.xlsx', sheet_name='SuppleTable2_v2.txt', skiprows=2)


# In[28]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[29]:


original_data.head()


# In[30]:


original_data['orf'] = original_data['ORF'].astype(str)


# In[31]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[32]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[33]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[34]:


original_data.set_index('orf', inplace=True)


# In[35]:


data_cols = [c for c in original_data.columns if 'yko' in c]
data_cols


# In[36]:


original_data = original_data[data_cols].copy()


# In[37]:


original_data = original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[38]:


original_data = original_data.groupby(original_data.index).mean()


# In[39]:


original_data.shape


# In[40]:


original_data.columns = [4711, 4781, 4782, 4783, 4711, 4781, 4782, 4783]


# In[41]:


original_data = original_data.T
original_data = original_data.groupby(original_data.index).mean()
original_data = original_data.T


# In[42]:


original_data.shape


# In[43]:


original_data.head()


# # Prepare the final dataset

# In[48]:


data = original_data.copy()


# In[49]:


dataset_ids = [4711, 4781, 4782, 4783]
datasets = datasets.reindex(index=dataset_ids)


# In[50]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[51]:


data.head()


# ## Subset to the genes currently in SGD

# In[52]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[53]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[54]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[55]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[56]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[57]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[58]:


from IO.save_data_to_db3 import *


# In[59]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




