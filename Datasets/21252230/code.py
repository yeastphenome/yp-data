#!/usr/bin/env python
# coding: utf-8

# In[2]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[3]:


paper_pmid = 21252230
paper_name = 'martin_cunningham_2011' 


# In[4]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[5]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[6]:


original_data = pd.read_excel('raw_data/1-s2.0-S0021925820538970-mmc1.xls', sheet_name='data')


# In[7]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[8]:


original_data.head()


# In[9]:


original_data['orf'] = original_data['ORF'].astype(str)


# In[10]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[11]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[13]:


original_data = original_data.loc[t,:]


# In[14]:


original_data.set_index('orf', inplace=True)


# In[15]:


data_cols = ['n(0)', 'n(FK)', 'n(TM)', 'n(TMFK)', 'n(AF)', 'n(AFFK)']


# In[16]:


original_data = original_data[data_cols].apply(pd.to_numeric, axis=1, errors='coerce')


# In[17]:


original_data = original_data.groupby(original_data.index).mean()


# In[18]:


original_data.shape


# In[20]:


# Normalize by control (n(0))
original_data = original_data.div(original_data['n(0)'], axis=0)


# In[22]:


original_data.drop(columns='n(0)', inplace=True)


# # Prepare the final dataset

# In[23]:


data = original_data.copy()


# In[24]:


dataset_ids = [21844, 21843, 21841, 21842, 21840]
datasets = datasets.reindex(index=dataset_ids)


# In[25]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[26]:


data.head()


# ## Subset to the genes currently in SGD

# In[27]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[28]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[29]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[30]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[31]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[32]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[33]:


from IO.save_data_to_db3 import *


# In[34]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




