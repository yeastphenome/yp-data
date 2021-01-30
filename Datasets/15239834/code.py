#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 15239834
paper_name = 'hartman_tippery_2004' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[52]:


original_data = pd.read_excel('raw_data/gb-2004-5-7-r49-s7.xlsx', sheet_name='data')


# In[53]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[54]:


original_data.head()


# In[55]:


original_data['orf'] = original_data['ORF'].astype(str)


# In[56]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[57]:


# Remove trailing "B" from ORFs
original_data['orf'] = original_data['orf'].apply(lambda x: x.strip('B') if '-' not in x else x)


# In[58]:


typo_fixes = {'YJR055WC':'YJR055W','YNL089CC':'YNL089C','YNL096CC':'YNL096C','YOR298C-AB':'YOR298C-A'}
original_data['orf'] = original_data['orf'].apply(lambda x: typo_fixes[x] if x in typo_fixes.keys() else x)


# In[59]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[60]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[61]:


original_data = original_data.loc[t,:]


# In[62]:


original_data.set_index('orf', inplace=True)


# In[63]:


original_data = original_data[['No HU- Growth Index','50mM HU Growth Index','150 mM HU Growth Index']].apply(pd.to_numeric, axis=1, errors='coerce')


# In[64]:


original_data.columns = ['unt','50','150']


# In[65]:


original_data['50'] = original_data['50'] - original_data['unt']
original_data['150'] = original_data['150'] - original_data['unt']


# In[66]:


original_data = original_data.groupby(original_data.index).mean()


# In[67]:


original_data.shape


# # Prepare the final dataset

# In[68]:


data = original_data.copy()


# In[69]:


dataset_ids = [16186, 52, 53]
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




