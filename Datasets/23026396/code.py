#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 23026396
paper_name = 'zhao_jiang_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[17]:


original_data = pd.read_excel('raw_data/mmc2.xlsx', sheet_name='Sheet1', skiprows=1)


# In[18]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[19]:


original_data.head()


# In[20]:


original_data['orf'] = original_data['&Systemic name'].astype(str)


# In[21]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[22]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[23]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[24]:


original_data = original_data.loc[t,:]


# In[25]:


original_data['data'] = -1


# In[26]:


original_data.set_index('orf', inplace=True)


# In[27]:


original_data = original_data[['data']].copy()


# In[28]:


original_data = original_data.groupby(original_data.index).mean()


# In[29]:


original_data.shape


# # Load & process tested strains

# In[32]:


tested = pd.read_excel('raw_data/DELETION LIBRARY.xlsx', sheet_name='DELETION LIBRARY', skiprows=1)


# In[33]:


tested.head()


# In[34]:


tested['orf'] = tested['ORF name'].astype(str)


# In[35]:


tested['orf'] = clean_orf(tested['orf'])


# In[38]:


tested.loc[tested['orf']=='YELOO1C','orf'] = 'YEL001C'


# In[39]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[40]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[41]:


tested = tested.loc[t,:]


# In[42]:


tested_orfs = tested['orf'].unique()


# In[43]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[44]:


tested_orfs = list(tested_orfs) + missing


# In[45]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[46]:


data = original_data.copy()


# In[47]:


dataset_ids = [1309]
datasets = datasets.reindex(index=dataset_ids)


# In[48]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[49]:


data.head()


# ## Subset to the genes currently in SGD

# In[50]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[51]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

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

data_all.head()


# # Print out

# In[55]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[56]:


from IO.save_data_to_db3 import *


# In[57]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




