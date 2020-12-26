#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 28404745
paper_name = 'dalton_conibear_2017' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/mc-E16-11-0772-s08.xlsx', sheet_name='master a alpha results')


# In[6]:


original_data.head()


# In[7]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[8]:


original_data['orf'] = original_data['ORF name'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[10]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[12]:


original_data.set_index('orf', inplace=True)


# In[13]:


for s in ['MAT α Z Score','MAT a Z Score']:
    original_data[s] = -pd.to_numeric(original_data[s], errors='coerce')


# In[14]:


original_data = original_data[['MAT α Z Score','MAT a Z Score']].copy()


# In[15]:


original_data = original_data.groupby(original_data.index).mean()


# In[16]:


original_data.shape


# In[17]:


dataset_ids = [16263, 16262]


# # Prepare the final dataset

# In[18]:


data = original_data.copy()


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


# In[24]:


data.head()


# # Normalize

# In[25]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[26]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[27]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[28]:


data_all.head()


# # Print out

# In[29]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[30]:


from IO.save_data_to_db3 import *


# In[31]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




