#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 25074250
paper_name = 'hwang_naganuma_2014' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[25]:


original_data1 = pd.read_excel('raw_data/tables_1_2.xlsx', sheet_name='Resistance')
original_data2 = pd.read_excel('raw_data/tables_1_2.xlsx', sheet_name='Sensitivity')


# In[26]:


original_data = pd.concat([original_data1, original_data2], axis=0)


# In[27]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[28]:


original_data.head()


# In[29]:


original_data['orfs'] = original_data['ORF'].astype(str)


# In[30]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[31]:


original_data = original_data.reset_index()


# In[32]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[33]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[34]:


# Subtracting 1, so that the non-hit data can be automatically assigned to 0 at future analysis steps
original_data['data'] = original_data['IC50 of deletion cells/IC50 of control cells'] - 1


# In[35]:


original_data.set_index('orfs', inplace=True)
original_data.index.name='orf'


# In[36]:


original_data = original_data[['data']].copy()


# In[37]:


original_data = original_data.groupby(original_data.index).mean()


# In[38]:


original_data.head()


# # Prepare the final dataset

# In[39]:


data = original_data.copy()


# In[40]:


dataset_ids = [16508]
datasets = datasets.reindex(index=dataset_ids)


# In[41]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[42]:


data.head()


# ## Subset to the genes currently in SGD

# In[43]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[44]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[46]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[47]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[48]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[49]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[50]:


from IO.save_data_to_db3 import *


# In[51]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




