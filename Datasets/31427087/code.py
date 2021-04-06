#!/usr/bin/env python
# coding: utf-8

# In[2]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[3]:


paper_pmid = 31427087
paper_name = 'zhao_han_2019' 


# In[4]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[5]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[75]:


original_data = pd.read_excel('raw_data/1-s2.0-S0006291X19314111-mmc3.xls', sheet_name='Sheet1', skiprows=1)


# In[76]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[77]:


original_data.head()


# In[78]:


original_data['orf'] = original_data['Systematic name'].astype(str)


# In[79]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[80]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[81]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[82]:


original_data['data'] = -1


# In[83]:


original_data.set_index('orf', inplace=True)


# In[84]:


original_data = original_data[['data']].copy()


# In[85]:


original_data = original_data.groupby(original_data.index).mean()


# In[86]:


original_data.shape


# In[46]:


original_data2 = pd.read_excel('raw_data/1-s2.0-S0006291X19314111-mmc4.xls', sheet_name='Sheet1', skiprows=1)


# In[47]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[48]:


original_data2.head()


# In[49]:


original_data2['orf'] = original_data2['Systematic name'].astype(str)


# In[50]:


# Eliminate all white spaces & capitalize
original_data2['orf'] = clean_orf(original_data2['orf'])


# In[51]:


# Translate to ORFs 
original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')


# In[52]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[53]:


original_data2['data'] = pd.to_numeric(original_data2['Resistant degree'], errors='coerce')


# In[54]:


original_data2.set_index('orf', inplace=True)


# In[55]:


original_data2 = original_data2[['data']].copy()


# In[56]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[57]:


original_data2.shape


# In[89]:


original_data3 = original_data.join(original_data2, how='outer', lsuffix='_sens', rsuffix='_res')


# In[90]:


original_data3[original_data3.isnull()] = 0


# In[92]:


original_data3.min()


# # Prepare the final dataset

# In[93]:


data = original_data3.copy()


# In[94]:


dataset_ids = [16709, 16711]
datasets = datasets.reindex(index=dataset_ids)


# In[95]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[96]:


data.head()


# ## Subset to the genes currently in SGD

# In[97]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[98]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[99]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[100]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[101]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[102]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[103]:


from IO.save_data_to_db3 import *


# In[104]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




