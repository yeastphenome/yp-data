#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 29783050
paper_name = 'jiang_papadopoulos_2018' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_csv('raw_data/hits.txt', sep='\t', header=None, names=['orfs'])


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data['orfs'] = original_data['orfs'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[10]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])


# In[12]:


print(original_data.loc[~t,])


# In[13]:


original_data['score'] = -1


# In[14]:


# If the same strain is present more than once, average its values
data = original_data.groupby('orfs')['score'].mean().to_frame()


# In[16]:


data.index.name='orf'


# In[15]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Prepare the final dataset

# In[18]:


dataset_ids = [16444]
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


# In[23]:


data.head()


# # Normalize

# In[24]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[25]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[26]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[27]:


data_all.head()


# # Print out

# In[28]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[29]:


from IO.save_data_to_db3 import *


# In[30]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




