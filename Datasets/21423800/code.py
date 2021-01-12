#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 21423800
paper_name = 'brett_rao_2011' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/journal.pone.0017619.s003.xlsx', sheet_name='Unsorted Data')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data['ORF name'] = original_data['ORF name'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['ORF name'] = clean_orf(original_data['ORF name'])


# In[9]:


# Translate to ORFs 
original_data['ORF name'] = translate_sc(original_data['ORF name'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF name'])
print(original_data.loc[~t,])


# In[11]:


original_data.set_index('ORF name', inplace=True)
original_data.index.name='orf'


# In[13]:


original_data = original_data[['Growth at pH 2.7','Growth at pH 4.0','Growth at pH 7.0','pHv at pH 2.7','pHv at pH 4.0','pHv at pH 7.0']].copy()


# In[14]:


original_data = original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[15]:


original_data = original_data.groupby(original_data.index).mean()


# In[16]:


original_data.shape


# # Prepare the final dataset

# In[17]:


data = original_data.copy()


# In[18]:


dataset_ids = [529,530,531,59,60,61]
datasets = datasets.reindex(index=dataset_ids)


# In[19]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[20]:


data.head()


# ## Subset to the genes currently in SGD

# In[21]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[22]:


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




