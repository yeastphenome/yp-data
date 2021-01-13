#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 20973990
paper_name = 'mira_sa_correia_2010' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[44]:


original_data = pd.read_excel('raw_data/1475-2859-9-79-s1.xlsx', sheet_name='Sheet1', skiprows=4)


# In[45]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[46]:


original_data.head()


# In[47]:


original_data['gene'] = original_data['Gene'].astype(str)


# In[48]:


# Eliminate all white spaces & capitalize
original_data['gene'] = clean_genename(original_data['gene'])


# In[49]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['gene'], to='orf')


# In[50]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[51]:


original_data = original_data.loc[t,:]


# In[52]:


original_data['data'] = original_data['Susceptibility to acetic acid2'].apply(lambda x: -len(x) if isinstance(x, str) else 0)


# In[53]:


original_data.set_index('orf', inplace=True)


# In[54]:


original_data = original_data[['data']].copy()


# In[55]:


original_data = original_data.groupby(original_data.index).mean()


# In[56]:


original_data.shape


# # Load & process tested strains

# In[57]:


tested = pd.read_excel('raw_data/List of strains tested.xlsx', sheet_name='Tabelle2')


# In[58]:


tested.head()


# In[59]:


tested['orf'] = tested['ORF'].astype(str)


# In[60]:


tested['orf'] = clean_orf(tested['orf'])


# In[61]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[62]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[63]:


tested_orfs = tested['orf'].unique()


# In[64]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[65]:


tested_orfs = list(tested_orfs) + missing


# In[66]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[67]:


data = original_data.copy()


# In[68]:


dataset_ids = [101]
datasets = datasets.reindex(index=dataset_ids)


# In[69]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[70]:


data.head()


# ## Subset to the genes currently in SGD

# In[71]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[72]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[73]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[74]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[75]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[76]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[77]:


from IO.save_data_to_db3 import *


# In[78]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




