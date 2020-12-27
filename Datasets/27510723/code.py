#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 27510723
paper_name = 'kwon_koo_2016' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data (1)

# In[31]:


original_data1 = pd.read_excel('raw_data/425_2016_2579_MOESM3_ESM.xlsx', sheet_name='Sheet1', skiprows=2)


# In[32]:


print('Original data dimensions: %d x %d' % (original_data1.shape))


# In[33]:


original_data1.head()


# In[34]:


original_data1['orf'] = original_data1['strain'].astype(str)


# In[35]:


# Eliminate all white spaces & capitalize
original_data1['orf'] = clean_orf(original_data1['orf'])


# In[36]:


# Translate to ORFs 
original_data1['orf'] = translate_sc(original_data1['orf'], to='orf')


# In[37]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orf'])
print(original_data1.loc[~t,])


# In[38]:


original_data1['data'] = -original_data1['LogRatio_BestTag_Zscore']


# In[39]:


original_data1.set_index('orf', inplace=True)


# In[40]:


original_data1 = original_data1[['data']].copy()


# In[41]:


original_data1 = original_data1.groupby(original_data1.index).mean()


# In[42]:


original_data1.shape


# # Load & process data (2)

# In[18]:


original_data2 = pd.read_excel('raw_data/het/Background Corrected Ratios And Robust Z-Scores (FD Scores)(all data).xlsx', sheet_name='Sheet1')


# In[20]:


original_data2['orf'] = original_data2['strain'].astype(str)


# In[21]:


# Eliminate all white spaces & capitalize
original_data2['orf'] = clean_orf(original_data2['orf'])


# In[22]:


# Translate to ORFs 
original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')


# In[23]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[26]:


original_data2['data'] = -pd.to_numeric(original_data2['LogRatio_BestTag_Zscore'], errors='coerce')


# In[27]:


original_data2.set_index('orf', inplace=True)


# In[28]:


original_data2 = original_data2[['data']].copy()


# In[29]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[30]:


original_data2.shape


# # Merge data

# In[43]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[44]:


original_data.head()


# # Prepare the final dataset

# In[45]:


data = original_data.copy()


# In[46]:


dataset_ids = [4826, 5181]
datasets = datasets.reindex(index=dataset_ids)


# In[47]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[48]:


data.head()


# ## Subset to the genes currently in SGD

# In[49]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[50]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[51]:


data.head()


# # Normalize

# In[52]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[53]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[54]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[55]:


data_all.head()


# # Print out

# In[56]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[57]:


from IO.save_data_to_db3 import *


# In[58]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




