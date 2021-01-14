#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 20206679
paper_name = 'zhao_jiang_2010' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[21]:


original_data = pd.read_excel('raw_data/hits.xlsx', sheet_name='Sheet1')


# In[22]:


original_data.columns = [x.strip() for x in original_data.columns]


# In[23]:


original_data.head()


# In[24]:


cols = ['Systemic Name','Systemic Name.1','Systemic Name.2','Systemic Name.3']
original_data = pd.concat([original_data[col] for col in cols], axis=0).to_frame()


# In[25]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[26]:


original_data.head()


# In[27]:


original_data['orfs'] = original_data[0].astype(str)


# In[28]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[29]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'].values, to='orf')


# In[30]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[31]:


original_data = original_data.loc[t,:]


# In[32]:


original_data.set_index('orfs', inplace=True)
original_data.index.name='orf'


# In[33]:


original_data['data'] = -1


# In[34]:


original_data = original_data[['data']].copy()


# In[35]:


original_data = original_data.groupby(original_data.index).mean()


# In[36]:


original_data.shape


# # Load & process tested strains

# In[38]:


tested = pd.read_excel('raw_data/DELETION LIBRARY.xlsx', sheet_name='DELETION LIBRARY', skiprows=1)


# In[39]:


tested.head()


# In[40]:


tested['orf'] = tested['ORF name'].astype(str)


# In[41]:


tested['orf'] = clean_orf(tested['orf'])


# In[42]:


tested.loc[tested['orf']=='YELOO1C','orf'] = 'YEL001C'


# In[43]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[44]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[45]:


tested = tested.loc[t,:]


# In[46]:


tested_orfs = tested['orf'].unique()


# In[47]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[48]:


tested_orfs = list(tested_orfs) + missing


# In[49]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# In[50]:


original_data.shape


# # Prepare the final dataset

# In[51]:


data = original_data.copy()


# In[52]:


dataset_ids = [16456]
datasets = datasets.reindex(index=dataset_ids)


# In[53]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[54]:


data.head()


# ## Subset to the genes currently in SGD

# In[55]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[56]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[58]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[59]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[60]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[62]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[63]:


from IO.save_data_to_db3 import *


# In[64]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




