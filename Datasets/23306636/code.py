#!/usr/bin/env python
# coding: utf-8

# In[2]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[3]:


paper_pmid = 23306636
paper_name = 'hoshida_akada_2013' 


# In[4]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[5]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[6]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, names=['genes'])


# In[7]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[8]:


original_data['genes'] = original_data['genes'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[10]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[12]:


original_data['data'] = 1


# In[13]:


original_data.set_index('orfs', inplace=True)
original_data.index.name='orf'


# In[14]:


original_data = original_data.groupby(original_data.index).mean()


# In[15]:


original_data.shape


# # Prepare the final dataset

# In[16]:


data = original_data[['data']].copy()


# In[17]:


dataset_ids = [16537]
datasets = datasets.reindex(index=dataset_ids)


# In[18]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[19]:


data.head()


# ## Subset to the genes currently in SGD

# In[20]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[21]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[22]:


data.head()


# # Normalize

# In[23]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[24]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[25]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[26]:


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



