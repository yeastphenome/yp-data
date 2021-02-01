#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 12972632
paper_name = 'huang_kolodner_2003' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[30]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, sep='\t')


# In[31]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[32]:


original_data.head()


# In[33]:


original_data['gene'] = original_data[0].astype(str)


# In[34]:


# Eliminate all white spaces & capitalize
original_data['gene'] = clean_genename(original_data['gene'])


# In[35]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['gene'], to='orf')


# In[36]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[37]:


original_data['data'] = pd.to_numeric(original_data[1], errors='coerce')


# In[38]:


original_data.set_index('orf', inplace=True)


# In[39]:


original_data['data'] = (original_data['data'] / original_data.loc['WILDTYPE','data']) - 1


# In[40]:


original_data = original_data[['data']].copy()


# In[41]:


original_data = original_data.groupby(original_data.index).mean()


# In[42]:


original_data.drop(index='WILDTYPE', inplace=True)


# In[43]:


original_data.shape


# # Prepare the final dataset

# In[44]:


data = original_data.copy()


# In[45]:


dataset_ids = [74]
datasets = datasets.reindex(index=dataset_ids)


# In[46]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[47]:


data.head()


# ## Subset to the genes currently in SGD

# In[48]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[49]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[50]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[51]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[52]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[53]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[54]:


from IO.save_data_to_db3 import *


# In[55]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




