#!/usr/bin/env python
# coding: utf-8

# In[43]:


get_ipython().run_line_magic('run', '../yp_utils.py')

import re


# # Initial setup

# In[2]:


paper_pmid = 15883361
paper_name = 'panavas_nagy_2005' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[30]:


original_data = pd.read_excel('raw_data/hits.xlsx', sheet_name='Sheet1')


# In[31]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[32]:


original_data.head()


# In[33]:


original_data['gene'] = original_data['Gene*'].astype(str)


# In[34]:


# Eliminate all white spaces & capitalize
original_data['gene'] = clean_genename(original_data['gene'])


# In[35]:


original_data['gene'] = original_data['gene'].apply(lambda x: x.replace(',',''))


# In[36]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['gene'], to='orf')


# In[37]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[38]:


original_data = original_data.loc[t,:]


# In[44]:


original_data['data'] = original_data['Replication†'].apply(lambda x: re.sub("[^0-9]", "", x))


# In[47]:


original_data['data'] = pd.to_numeric(original_data['data'], errors='coerce')


# In[48]:


original_data.set_index('orf', inplace=True)


# In[49]:


original_data = original_data[['data']].copy()


# In[50]:


original_data = original_data.groupby(original_data.index).mean()


# In[51]:


original_data.shape


# In[52]:


original_data.head()


# # Prepare the final dataset

# In[53]:


data = original_data.copy()


# In[54]:


dataset_ids = [16007]
datasets = datasets.reindex(index=dataset_ids)


# In[55]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[56]:


data.head()


# ## Subset to the genes currently in SGD

# In[57]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[58]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[59]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[60]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[61]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[62]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[63]:


from IO.save_data_to_db3 import *


# In[64]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




