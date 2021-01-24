#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 17308930
paper_name = 'wang_zhou_2007' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[53]:


original_data = pd.read_excel('raw_data/mutants sensitive to metal scarcity.xlsx', sheet_name='Sheet1')


# In[54]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[55]:


original_data.head()


# In[56]:


original_data['orf'] = original_data['ORF'].astype(str)


# In[57]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[58]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[59]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[60]:


data_switch = {'S': -1, 'S, M': -1, 'HS': -2, 'HS, M': -2}
original_data['data'] = original_data['Sensitivity to EDTA'].apply(lambda x: data_switch[x])


# In[61]:


original_data.set_index('orf', inplace=True)


# In[62]:


original_data = original_data[['data']].copy()


# In[63]:


original_data = original_data.groupby(original_data.index).mean()


# In[64]:


original_data.shape


# # Load data (2)

# In[65]:


original_data2 = pd.read_excel('raw_data/mutants sensitive to metal excess.xlsx', sheet_name='Sheet1')
original_data2.head()


# In[66]:


original_data2['orf'] = original_data2['ORF'].astype(str)


# In[67]:


original_data2['orf'] = clean_orf(original_data2['orf'])


# In[68]:


original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')


# In[69]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[70]:


original_data2.set_index('orf', inplace=True)


# In[71]:


original_data2 = original_data2[['Senitivity to Cu','Senitivity to Fe','Senitivity to Mn','Senitivity to Zn']].copy()


# In[72]:


data_switch = {'－': 0, 'S': -1, 'S, M': -1, 'HS': -2, 'HS, M': -2, 'HS,M': -2}
for c in original_data2.columns:
    original_data2[c] = original_data2[c].apply(lambda x: data_switch[x])


# In[73]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[74]:


original_data2.shape


# In[75]:


original_data2.head()


# # Merge

# In[76]:


original_data = original_data.join(original_data2, how='outer')


# In[77]:


original_data.shape


# In[78]:


original_data[original_data.isnull()] = 0


# In[79]:


original_data.head()


# # Prepare the final dataset

# In[80]:


data = original_data.copy()


# In[81]:


dataset_ids = [5246, 5247, 5248, 5249, 5250]
datasets = datasets.reindex(index=dataset_ids)


# In[82]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[83]:


data.head()


# ## Subset to the genes currently in SGD

# In[84]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[85]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[86]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[87]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[88]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

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




