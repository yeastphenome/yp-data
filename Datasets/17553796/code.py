#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 17553796
paper_name = 'thevissen_francois_2007' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[45]:


original_data = pd.read_csv('raw_data/hits_orfs_scores.txt', sep='\t', header=None)


# In[46]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[47]:


original_data.head()


# In[48]:


original_data['orf'] = original_data[0].astype(str)


# In[49]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[50]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[51]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[52]:


original_data['data'] = -original_data[1].astype(int)


# In[53]:


original_data.set_index('orf', inplace=True)


# In[54]:


original_data = original_data[['data']].copy()


# In[55]:


original_data = original_data.groupby(original_data.index).mean()


# In[56]:


original_data.shape


# # Load & process tested strains

# In[57]:


tested = pd.read_excel('raw_data/Euroscarf library.xlsx', sheet_name='Tabelle1')


# In[58]:


tested.head()


# In[59]:


tested['orf'] = tested['ORF'].astype(str)


# In[60]:


tested['orf'] = clean_orf(tested['orf'])


# In[61]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[62]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[63]:


tested = tested.loc[t,:]


# In[64]:


tested_orfs = tested['orf'].unique()


# In[65]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[66]:


tested_orfs = list(tested_orfs) + missing


# In[67]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[68]:


data = original_data.copy()


# In[69]:


dataset_ids = [116]
datasets = datasets.reindex(index=dataset_ids)


# In[70]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[71]:


data.head()


# ## Subset to the genes currently in SGD

# In[72]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[73]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[74]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[75]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[76]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[77]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[78]:


from IO.save_data_to_db3 import *


# In[79]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




