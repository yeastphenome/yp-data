#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 25773006
paper_name = 'du_jiang_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[34]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, names=['genes'])


# In[35]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[36]:


original_data['genes'] = original_data['genes'].astype(str)


# In[37]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[38]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[39]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])


# In[40]:


print(original_data.loc[~t,])


# In[41]:


original_data['data'] = -1


# In[42]:


original_data.set_index('orfs', inplace=True)
original_data.index.name = 'orf'


# In[43]:


original_data = original_data[['data']].copy()


# In[44]:


original_data = original_data.groupby(original_data.index).mean()


# In[45]:


original_data.shape


# # Load & process tested strains

# In[46]:


tested = pd.read_excel('raw_data/DELETION LIBRARY.xlsx', sheet_name='DELETION LIBRARY', skiprows=1)
tested['orf'] = tested['ORF name'].astype(str)
tested['orf'] = clean_orf(tested['orf'])
tested.loc[tested['orf'] == 'YELOO1C','orf'] = 'YEL001C'
tested['orf'] = translate_sc(tested['orf'], to='orf')
# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t])


# In[47]:


tested = tested.loc[t,:]


# In[48]:


tested_orfs = np.unique(tested['orf'].values)


# In[49]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[50]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[51]:


data = original_data.copy()


# In[52]:


dataset_ids = [16461]
datasets = datasets.reindex(index=dataset_ids)


# In[53]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[54]:


data.head()


# ## Subset to the genes currently in SGD

# In[55]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[56]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[57]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[58]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[59]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[60]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[61]:


from IO.save_data_to_db3 import *


# In[62]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




