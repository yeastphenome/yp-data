#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 28076367
paper_name = 'bozaquel_morais_montero_lomeli_2017' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[50]:


original_data = pd.read_excel('raw_data/journal.pone.0169682.s002.xlsx', sheet_name='PRIMARY HITS', skiprows=1)


# In[51]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[52]:


original_data.head()


# In[53]:


original_data['orf'] = original_data['SYSTEMATIC'].astype(str)


# In[54]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[55]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[56]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[57]:


original_data = original_data.loc[t,:]


# In[58]:


original_data['data1'] = -original_data.iloc[:,2:4].sum(axis=1)
original_data['data2'] = -original_data.iloc[:,4:6].sum(axis=1)


# In[59]:


original_data.set_index('orf', inplace=True)


# In[60]:


original_data = original_data[['data1','data2']].copy()


# In[61]:


original_data = original_data.groupby(original_data.index).mean()


# In[62]:


original_data.shape


# # Load & process tested strains

# In[63]:


tested = pd.read_excel('raw_data/Mat a collection YSC1053 open biosystems.xlsx', sheet_name='DATA')


# In[64]:


tested.head()


# In[65]:


tested['orf'] = tested['ORF name'].astype(str)


# In[66]:


tested['orf'] = clean_orf(tested['orf'])


# In[67]:


tested.loc[tested['orf']=='YLR287-A','orf'] = 'YLR287C-A'


# In[68]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[69]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[70]:


tested = tested.loc[t,:]


# In[71]:


tested_orfs = np.unique(tested['orf'].values)


# In[73]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[74]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# In[75]:


original_data.head()


# # Prepare the final dataset

# In[76]:


data = original_data.copy()


# In[77]:


dataset_ids = [16260, 16261]
datasets = datasets.reindex(index=dataset_ids)


# In[78]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[79]:


data.head()


# ## Subset to the genes currently in SGD

# In[80]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[81]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[82]:


data.head()


# # Normalize

# In[83]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[84]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[85]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[86]:


data_all.head()


# # Print out

# In[87]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[88]:


from IO.save_data_to_db3 import *


# In[89]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




