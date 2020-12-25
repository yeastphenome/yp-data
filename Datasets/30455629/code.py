#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 30455629
paper_name = 'fruhmann_cullin_2018' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/Table_1_The Impact of ESCRT on Aβ1-42 Induced Membrane Lesions in a Yeast Model for Alzheimer’s Disease.XLS', 
                              sheet_name='Feuil1', skiprows=2)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data['orf'] = original_data['systematic name'].astype(str)


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


original_data['data'] = 0


# In[13]:


# Enhancers of AB-42 toxicity = more severe growth defect
ix = original_data['Growth'] <= 2
original_data.loc[ix,'data'] = original_data.loc[ix,'Growth'] - 3

# Suppressors of AB-42 toxicity = less severe growth defect
ix = original_data['Growth'] >= 3
original_data.loc[ix,'data'] = original_data.loc[ix,'Growth'] - 2


# In[14]:


original_data.set_index('orf', inplace=True)


# In[15]:


original_data = original_data.groupby(original_data.index).mean()


# In[16]:


original_data.shape


# # Prepare the final dataset

# In[18]:


data = original_data[['data']].copy()


# In[19]:


dataset_ids = [16614]
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


data_norm = normalize_phenotypic_scores(data, has_tested=False)


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

# In[36]:


from IO.save_data_to_db2 import *


# In[37]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[38]:


save_data_to_db(data, paper_pmid)


# In[ ]:




