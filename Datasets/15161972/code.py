#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 15161972
paper_name = 'askree_mceachern_2004' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[41]:


original_data = pd.read_excel('raw_data/01263Table3.xlsx', sheet_name='Sheet1', skiprows=1)


# In[42]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[43]:


original_data.head()


# In[44]:


original_data['orf'] = original_data['ORF name'].astype(str)


# In[45]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[46]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[47]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[48]:


original_data = original_data.loc[t,:]


# In[49]:


original_data.set_index('orf', inplace=True)


# In[50]:


original_data = original_data[['Measurement 1','Measurement 2']].apply(pd.to_numeric, axis=1, errors='coerce').mean(axis=1).to_frame()


# In[51]:


original_data = original_data.groupby(original_data.index).mean()


# In[52]:


original_data.shape


# # Load & process tested strains

# In[53]:


tested = pd.read_excel('raw_data/S.cerftpmata.xlsx', sheet_name='Sheet1', skiprows=2)


# In[54]:


tested.head()


# In[55]:


tested['orf'] = tested['ORF name'].astype(str)


# In[56]:


tested['orf'] = clean_orf(tested['orf'])


# In[57]:


typo_fixes = {'YOLO57W':'YOL057W','YOLO62C':'YOL062C','YKLO72W':'YKL072W','YJL206-A':'YJL206C-A','YLR287-A':'YLR287C-A'}
tested['orf'] = tested['orf'].apply(lambda x: typo_fixes[x] if x in typo_fixes.keys() else x)


# In[58]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[59]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[60]:


tested = tested.loc[t,:]


# In[61]:


tested_orfs = tested['orf'].unique()


# In[62]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[63]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# In[64]:


original_data.head()


# # Prepare the final dataset

# In[65]:


data = original_data.copy()


# In[66]:


dataset_ids = [113]
datasets = datasets.reindex(index=dataset_ids)


# In[67]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[68]:


data.head()


# ## Subset to the genes currently in SGD

# In[69]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[70]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[71]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[72]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[73]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[74]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[75]:


from IO.save_data_to_db3 import *


# In[76]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




