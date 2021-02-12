#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 22384326
paper_name = 'kloimwieder_winston_2011' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[40]:


original_data = pd.read_excel('raw_data/TableS1.xlsx', sheet_name='Table 1', skiprows=2)


# In[41]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[42]:


original_data.head()


# In[43]:


orfs = pd.concat([original_data.iloc[:,1], 
                  original_data.iloc[:,5],
                  original_data.iloc[:,9]], axis=0)


# In[44]:


original_data = pd.DataFrame(data={'orf': orfs})


# In[45]:


original_data['orf'] = original_data['orf'].astype(str)


# In[46]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[47]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'].values, to='orf')


# In[48]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[50]:


original_data = original_data.loc[t,:]


# In[51]:


original_data['data'] = -1


# In[52]:


original_data.set_index('orf', inplace=True)


# In[53]:


original_data = original_data[['data']].copy()


# In[54]:


original_data = original_data.groupby(original_data.index).mean()


# In[55]:


original_data.shape


# # Load & process tested strains

# In[57]:


tested = pd.read_excel('raw_data/Homo_diploids_101501.xlsx', sheet_name='Homo_diploids_101501.txt', skiprows=1)


# In[58]:


tested.head()


# In[59]:


tested['orf'] = tested['ORF name'].astype(str)


# In[60]:


tested['orf'] = clean_orf(tested['orf'])


# In[61]:


tested.loc[tested['orf']=='YELOO1C','orf'] = 'YEL001C'


# In[62]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[63]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[64]:


tested = tested.loc[t,:]


# In[65]:


tested_orfs = tested['orf'].unique()


# In[66]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[67]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[68]:


data = original_data.copy()


# In[69]:


dataset_ids = [16136]
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




