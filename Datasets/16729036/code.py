#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 16729036
paper_name = 'dudley_church_2005' 


# In[46]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[47]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/pheno_data.xlsx', sheet_name='YPDNorm')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data['orf'] = original_data['YName'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[14]:


original_data.loc[original_data['orf']=='YELOO1C','orf'] = 'YEL001C'


# In[15]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[16]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[17]:


original_data = original_data.loc[t,:]


# In[18]:


original_data.set_index('orf', inplace=True)


# In[25]:


data_cols = original_data.columns[9:53].values


# In[26]:


original_data = original_data[data_cols].apply(pd.to_numeric, axis=1, errors='coerce')


# In[27]:


original_data = original_data.groupby(original_data.index).mean()


# In[28]:


original_data.shape


# # Load & process dataset ids

# In[32]:


dt = pd.read_excel('extras/drug_dataset.xlsx', sheet_name='Sheet1', header=None)


# In[33]:


dt.head()


# In[34]:


dt.set_index(0, inplace=True)


# In[35]:


dt = dt.reindex(index=original_data.columns.values)


# In[37]:


dataset_ids = dt[1].values


# In[42]:


original_data.columns = dataset_ids


# In[43]:


original_data = original_data.T
original_data = original_data.groupby(original_data.index).mean()
original_data = original_data.T


# # Prepare the final dataset

# In[48]:


data = original_data.copy()


# In[49]:


dataset_ids = original_data.columns.values
datasets = datasets.reindex(index=dataset_ids)


# In[50]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[51]:


data.head()


# ## Subset to the genes currently in SGD

# In[52]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[53]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[54]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[55]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[56]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[57]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[58]:


from IO.save_data_to_db3 import *


# In[59]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




