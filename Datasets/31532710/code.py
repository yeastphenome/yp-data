#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 31532710
paper_name = 'goke_walter_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/Tables_S1-1_S1-2.xlsx', sheet_name='Table 1', skiprows=1)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[22]:


orfs1 = pd.concat((original_data.iloc[:,0], original_data.iloc[:,1], original_data.iloc[:,2]), axis=0).unique()


# In[23]:


orfs2 = pd.concat((original_data.iloc[:,4], original_data.iloc[:,5], original_data.iloc[:,6], original_data.iloc[:,7]), axis=0).unique()


# In[24]:


orfs3 = pd.concat((original_data.iloc[:,9], original_data.iloc[:,10]), axis=0).unique()


# In[25]:


orfs4 = pd.concat((original_data.iloc[:,12], original_data.iloc[:,13], 
                   original_data.iloc[:,14], original_data.iloc[:,15],original_data.iloc[:,16]), axis=0).unique()


# In[26]:


orfs1 = pd.DataFrame(index=orfs1, data={'data': 1})
orfs2 = pd.DataFrame(index=orfs2, data={'data': 1})
orfs3 = pd.DataFrame(index=orfs3, data={'data': 1})
orfs4 = pd.DataFrame(index=orfs4, data={'data': -1})


# In[27]:


original_data = pd.concat([orfs1, orfs2, orfs3, orfs4], axis=1)


# In[29]:


original_data[original_data.isnull()] = 0


# In[30]:


original_data['data_all'] = original_data.sum(axis=1)


# In[34]:


original_data = original_data.reset_index()


# In[36]:


original_data['orf'] = original_data['index'].astype(str)


# In[37]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[38]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[39]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[40]:


original_data = original_data.loc[t,:]


# In[41]:


original_data.set_index('orf', inplace=True)


# In[42]:


original_data = original_data[['data_all']].copy()


# In[43]:


original_data = original_data.groupby(original_data.index).mean()


# In[44]:


original_data.shape


# # Prepare the final dataset

# In[45]:


data = original_data.copy()


# In[46]:


dataset_ids = [21878]
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

data.head()


# # Normalize

# In[51]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[52]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[53]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[54]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[55]:


from IO.save_data_to_db3 import *


# In[56]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




