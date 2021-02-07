#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 11907261
paper_name = 'seeley_eitzen_2002' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[15]:


original_data = pd.read_excel('raw_data/data.xlsx', sheet_name='Sheet1')


# In[16]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[17]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[18]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[19]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'], to='orf')


# In[20]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[21]:


original_data.set_index('ORF', inplace=True)
original_data.index.name = 'orf'


# In[22]:


original_data = original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[23]:


original_data[original_data.isnull()] = 0


# In[24]:


original_data = original_data.groupby(original_data.index).mean()


# In[25]:


original_data.shape


# # Prepare the final dataset

# In[26]:


data = original_data.copy()


# In[27]:


dataset_ids = [16542, 16539, 16540, 16541]
datasets = datasets.reindex(index=dataset_ids)


# In[28]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[29]:


data.head()


# ## Subset to the genes currently in SGD

# In[30]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[31]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[32]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[33]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[34]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[35]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[36]:


from IO.save_data_to_db3 import *


# In[38]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




