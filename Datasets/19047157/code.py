#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 19047157
paper_name = 'herrero_moreno_2008' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[23]:


original_data1 = pd.read_excel('raw_data/s_table_1.xlsx', sheet_name='Table 1', skiprows=2)


# In[24]:


print('Original data dimensions: %d x %d' % (original_data1.shape))


# In[25]:


original_data1.head()


# In[26]:


orf_cols = ['Unnamed: 1','Unnamed: 4','Unnamed: 7', 'Unnamed: 10']
lst = [original_data1[c] for c in orf_cols]
original_data1 = pd.concat(lst, axis=0, ignore_index=True)
original_data1.shape


# In[28]:


original_data1 = original_data1.to_frame()
original_data1.head()


# In[29]:


original_data1['orf'] = original_data1[0].astype(str)


# In[30]:


# Eliminate all white spaces & capitalize
original_data1['orf'] = clean_orf(original_data1['orf'])


# In[31]:


# Translate to ORFs 
original_data1['orf'] = translate_sc(original_data1['orf'], to='orf')


# In[32]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orf'])
print(original_data1.loc[~t,])


# In[33]:


original_data1 = original_data1.loc[t,:]


# In[34]:


original_data1['data'] = -1


# In[36]:


original_data1.set_index('orf', inplace=True)


# In[37]:


original_data1 = original_data1[['data']].copy()


# In[38]:


original_data1 = original_data1.groupby(original_data1.index).mean()


# In[39]:


original_data1.shape


# # Load data (2)

# In[40]:


original_data2 = pd.read_excel('raw_data/stab_3.xlsx', sheet_name='Table 1', skiprows=2)


# In[41]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[42]:


original_data2.head()


# In[43]:


original_data2['orf'] = original_data2['ORF'].astype(str)


# In[44]:


original_data2['orf'] = clean_orf(original_data2['orf'])


# In[45]:


original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')


# In[46]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[47]:


original_data2['data'] = 1


# In[48]:


original_data2.set_index('orf', inplace=True)


# In[49]:


original_data2 = original_data2[['data']].copy()


# In[50]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[51]:


original_data2.shape


# # Merge

# In[52]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[54]:


original_data['data'] = original_data.mean(axis=1)


# In[57]:


original_data.drop(columns=['data_1','data_2'], inplace=True)


# # Prepare the final dataset

# In[58]:


data = original_data.copy()


# In[59]:


dataset_ids = [106]
datasets = datasets.reindex(index=dataset_ids)


# In[60]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[61]:


data.head()


# ## Subset to the genes currently in SGD

# In[62]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[63]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[64]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[65]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[66]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[67]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[68]:


from IO.save_data_to_db3 import *


# In[69]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




