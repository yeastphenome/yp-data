#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 26646153
paper_name = 'choy_basrai_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[26]:


original_data_list = []
for p in np.arange(14)+1:
    dt = pd.read_excel('raw_data/supp_g3.115.022244_TableS2.xlsx', sheet_name='S1 P' + str(p))
    dt.columns = [c.lower() for c in dt.columns]
    original_data_list.append(dt)


# In[27]:


for p in np.arange(14)+1:
    dt = pd.read_excel('raw_data/supp_g3.115.022244_TableS3.xlsx', sheet_name='S2 P' + str(p))
    dt.columns = [c.lower() for c in dt.columns]
    original_data_list.append(dt)


# In[28]:


len(original_data_list)


# In[29]:


original_data = pd.concat(original_data_list, axis=0, ignore_index=True)


# In[30]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[31]:


original_data.head()


# In[32]:


original_data['orf'] = original_data['array orf'].astype(str)


# In[33]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[34]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[35]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[37]:


original_data = original_data.loc[t,:]


# In[38]:


original_data['data'] = original_data['ratio']


# In[39]:


original_data.set_index('orf', inplace=True)


# In[40]:


original_data = original_data[['data']].copy()


# In[41]:


original_data = original_data.groupby(original_data.index).mean()


# In[42]:


original_data.shape


# # Prepare the final dataset

# In[43]:


data = original_data.copy()


# In[44]:


dataset_ids = [768]
datasets = datasets.reindex(index=dataset_ids)


# In[45]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[46]:


data.head()


# ## Subset to the genes currently in SGD

# In[47]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[48]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[49]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[50]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[51]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[52]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[53]:


from IO.save_data_to_db3 import *


# In[54]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:



