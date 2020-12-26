#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 28209762
paper_name = 'bae_swanson_2017' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[52]:


original_data = pd.read_excel('raw_data/Table1.xlsx', sheet_name='Sheet1')


# In[53]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[54]:


original_data.head()


# In[55]:


original_data['orf'] = original_data['ORF ID'].astype(str)


# In[56]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[57]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[58]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[59]:


original_data = original_data.loc[t,:]


# In[60]:


original_data['data'] = -1


# In[61]:


original_data.set_index('orf', inplace=True)


# In[62]:


original_data = original_data.groupby(original_data.index).mean()


# In[63]:


original_data.shape


# In[64]:


dataset_ids = [11809]


# # Load & process tested strains

# In[65]:


tested = pd.read_excel('raw_data/Heterozygous_diploid_obs_v7.0.xlsx', sheet_name='Heterozygous Diploid_obs')


# In[66]:


tested.head()


# In[67]:


tested['orf'] = tested['ORF name'].astype(str)


# In[68]:


# Remove the underscore annotations
tested['orf'] = tested['orf'].apply(lambda x: x.split('_')[0])


# In[69]:


tested['orf'] = clean_orf(tested['orf'])


# In[70]:


typo_fixes = {'YCLO51W':'YCL051W','YHR139C-':'YHR139C-A','YGR122C-':'YGR122C-A'}
for s in typo_fixes.keys():
    tested.loc[tested['orf']==s,'orf'] = typo_fixes[s]


# In[71]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[72]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[73]:


tested_orfs = np.unique(tested['orf'].values)


# In[74]:


tested_orfs.shape


# In[75]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[76]:


data = original_data[['data']].copy()


# In[77]:


dataset_ids = [11809]
datasets = datasets.reindex(index=dataset_ids)


# In[78]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[79]:


data.head()


# In[80]:


data.sum()


# ## Subset to the genes currently in SGD

# In[81]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[82]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[83]:


data.head()


# # Normalize

# In[84]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[85]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[86]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[87]:


data_all.head()


# # Print out

# In[88]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[89]:


from IO.save_data_to_db3 import *


# In[90]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




