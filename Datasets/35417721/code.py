#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 35417721
paper_name = 'stolp_hardwick_2022' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[41]:


original_data = pd.read_excel('raw_data/1-s2.0-S2211124722003990-mmc2.xlsx', sheet_name='TableS1_Heat-ramp screen', 
                              skiprows=1, skipfooter=25)


# In[42]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[43]:


original_data.head()


# In[44]:


original_data['orf'] = original_data['orf_name'].astype(str)


# In[45]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[46]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[47]:


original_data.loc[original_data['orf']=='YLR287-A','orf'] = 'YLR287C-A'


# In[48]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[49]:


original_data['data'] = original_data['corr_z']


# In[50]:


original_data.set_index('orf', inplace=True)


# In[51]:


original_data = original_data[['data']].copy()


# In[52]:


original_data = original_data.groupby(original_data.index).mean()


# In[53]:


original_data.shape


# # Prepare the final dataset

# In[54]:


data = original_data.copy()


# In[55]:


dataset_ids = [22220]
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



