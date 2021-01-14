#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 20337531
paper_name = 'dias_sa_correia_2010' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[43]:


original_data1 = pd.read_csv('raw_data/TableS1.txt', header=None, names=['genes','orfs'], sep='\t')
original_data2 = pd.read_csv('raw_data/TableS2.txt', header=None, names=['genes','orfs'], sep='\t')


# In[44]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[45]:


original_data1['orfs'] = original_data1['orfs'].astype(str)
original_data2['orfs'] = original_data2['orfs'].astype(str)


# In[46]:


# Eliminate all white spaces & capitalize
original_data1['orfs'] = clean_orf(original_data1['orfs'])
original_data2['orfs'] = clean_orf(original_data2['orfs'])


# In[47]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['orfs'], to='orf')
original_data2['orfs'] = translate_sc(original_data2['orfs'], to='orf')


# In[48]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[49]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[50]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)
original_data1.index.name='orf'
original_data2.index.name='orf'


# In[51]:


# High-susceptibility
original_data1['data'] = -2

# Moderate-susceptibility
original_data2['data'] = -1


# In[52]:


original_data1 = original_data1[['data']].copy()
original_data2 = original_data2[['data']].copy()


# In[53]:


original_data1 = original_data1.groupby(original_data1.index).mean()
original_data2 = original_data2.groupby(original_data2.index).mean()


# In[54]:


original_data1.shape


# In[55]:


original_data2.shape


# In[56]:


original_data = original_data1.join(original_data2, lsuffix='_high', rsuffix='_mod', how='outer')


# In[57]:


original_data['data'] = original_data[['data_high','data_mod']].mean(axis=1)


# In[58]:


original_data.drop(columns=['data_high','data_mod'], inplace=True)


# In[59]:


original_data.shape


# # Load & process tested strains

# In[60]:


tested = pd.read_excel('raw_data/List of strains tested.xlsx', sheet_name='Tabelle2')


# In[61]:


tested.head()


# In[62]:


tested['orf'] = tested['ORF'].astype(str)


# In[63]:


tested['orf'] = clean_orf(tested['orf'])


# In[64]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[65]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,:])


# In[66]:


tested_orfs = tested['orf'].unique()


# In[67]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[68]:


# Remove missing (the data list contains HAP and HET screen results, so some of the genes are likely to be essential)
original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# In[69]:


original_data.shape


# # Prepare the final dataset

# In[70]:


data = original_data.copy()


# In[71]:


dataset_ids = [16607]
datasets = datasets.reindex(index=dataset_ids)


# In[73]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[74]:


data.head()


# ## Subset to the genes currently in SGD

# In[75]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[76]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[77]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[78]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[79]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[80]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[81]:


from IO.save_data_to_db3 import *


# In[82]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




