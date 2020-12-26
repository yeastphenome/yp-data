#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 28428821
paper_name = 'henriques_sa_correia_2017' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data (part 1)

# In[40]:


original_data1 = pd.read_excel('raw_data/13068_2017_781_MOESM1_ESM.xlsx', sheet_name='Table S1', skiprows=2)
original_data1.head()


# In[41]:


print('Original data dimensions: %d x %d' % (original_data1.shape))


# In[42]:


original_data1['genes'] = original_data1.iloc[:,0].astype(str)


# In[43]:


# Eliminate all white spaces & capitalize
original_data1['genes'] = clean_genename(original_data1['genes'])


# In[44]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['genes'], to='orf')


# In[45]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[46]:


original_data1 = original_data1.loc[t,:]


# In[47]:


original_data1['data'] = original_data1.iloc[:,2].apply(lambda x: -len(x))


# In[48]:


original_data1.set_index('orfs', inplace=True)
original_data1.index.name = 'orf'


# In[49]:


original_data1 = original_data1[['data']].copy()


# In[50]:


original_data1 = original_data1.groupby(original_data1.index).mean()


# In[51]:


original_data1.shape


# # Load & process the data (part 2)

# In[52]:


original_data2 = pd.read_excel('raw_data/13068_2017_781_MOESM2_ESM.xlsx', sheet_name='Table S2', skiprows=2)
original_data2.head()


# In[54]:


original_data2['genes'] = original_data2.iloc[:,0].astype(str)


# In[55]:


# Eliminate all white spaces & capitalize
original_data2['genes'] = clean_genename(original_data2['genes'])


# In[56]:


# Translate to ORFs 
original_data2['orfs'] = translate_sc(original_data2['genes'], to='orf')


# In[57]:


## Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[58]:


original_data2 = original_data2.loc[t,:]


# In[59]:


original_data2['data'] = 1


# In[60]:


original_data2.set_index('orfs', inplace=True)
original_data2.index.name = 'orf'


# In[61]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[62]:


original_data2.shape


# # Merge

# In[63]:


original_data = pd.concat((original_data1, original_data2), axis=0)


# In[66]:


original_data = original_data.groupby(original_data.index).mean()


# In[67]:


original_data.shape


# # Load & process tested strains

# In[68]:


tested = pd.read_excel('raw_data/List of strains tested.xlsx', sheet_name='Tabelle2')


# In[69]:


tested['orf'] = tested['ORF'].astype(str)


# In[70]:


tested['orf'] = clean_orf(tested['orf'])


# In[71]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[72]:


t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[73]:


tested_orfs = np.unique(tested['orf'].values)


# In[74]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[75]:


tested_orfs = list(tested_orfs) + missing


# In[76]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# In[77]:


original_data.shape


# # Prepare the final dataset

# In[78]:


data = original_data[['data']].copy()


# In[79]:


dataset_ids = [16264]
datasets = datasets.reindex(index=dataset_ids)


# In[80]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[81]:


data.head()


# ## Subset to the genes currently in SGD

# In[82]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[83]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[84]:


data.head()


# # Normalize

# In[85]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[86]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[87]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[88]:


data_all.head()


# # Print out

# In[89]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[90]:


from IO.save_data_to_db3 import *


# In[91]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




