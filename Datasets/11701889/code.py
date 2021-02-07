#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[47]:


paper_pmid = 11701889
paper_name = 'ooi_boeke_2001' 


# In[48]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[49]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[50]:


original_data = pd.read_excel('raw_data/ooi_boeke_2001.xlsx', sheet_name='Sheet1', header=None)


# In[51]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[52]:


original_data.head()


# In[53]:


original_data['gene'] = original_data[0].astype(str)


# In[54]:


# Eliminate all white spaces & capitalize
original_data['gene'] = clean_genename(original_data['gene'])


# In[55]:


original_data.loc[original_data['gene']=='GPE2','gene'] = 'GPB2'


# In[56]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['gene'], to='orf')


# In[57]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[58]:


data_switch = {'-': 0, '+': -1}
for c in [1,2]:
    original_data[c] = original_data[c].apply(lambda x: data_switch[x])


# In[59]:


original_data.set_index('orf', inplace=True)


# In[60]:


original_data = original_data[[1,2]].copy()


# In[61]:


original_data = original_data.groupby(original_data.index).mean()


# In[62]:


original_data.shape


# # Prepare the final dataset

# In[63]:


data = original_data.copy()


# In[64]:


dataset_ids = [315, 2]
datasets = datasets.reindex(index=dataset_ids)


# In[65]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[66]:


data.head()


# ## Subset to the genes currently in SGD

# In[67]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[68]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[69]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[70]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[71]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[72]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[73]:


from IO.save_data_to_db3 import *


# In[74]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




