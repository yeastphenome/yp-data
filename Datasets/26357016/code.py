#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 26357016
paper_name = 'frohlich_walther_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/Myriocin screen.xlsx', sheet_name='Myriocin screen')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


# Keep only deletions
original_data = original_data.loc[original_data['Mutation']=='DELETION',:]


# In[8]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[10]:


original_data.loc[original_data['ORF']=='YLR287-A','ORF'] = 'YLR287C-A'


# In[11]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'].values, to='orf')


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[13]:


# Normalize by untreated
original_data = original_data[['ORF','WT.mean','MYR.mean']].copy()
original_data['MYR.mean'] = original_data['MYR.mean'] / original_data['WT.mean']


# In[14]:


original_data.set_index('ORF', inplace=True)
original_data.index.name='orf'


# In[15]:


original_data = original_data.groupby(original_data.index).mean()


# In[16]:


original_data.shape


# # Prepare the final dataset

# In[17]:


data = original_data.copy()


# In[18]:


dataset_ids = [16496,16458]


# In[19]:


datasets = datasets.reindex(index=dataset_ids)


# In[20]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[21]:


data.head()


# ## Subset to the genes currently in SGD

# In[22]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[23]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[24]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[25]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[26]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[27]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[28]:


from IO.save_data_to_db3 import *


# In[29]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




