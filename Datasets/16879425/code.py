#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 16879425
paper_name = 'fujita_iwahashi_2006' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


files = ['Table1.txt','Table2.txt','Table3.txt']


# In[11]:


original_data_list = []
for f in files:
    original_data = pd.read_csv('raw_data/' + f, sep='\t', skiprows=1)
    print('Original data dimensions: %d x %d' % (original_data.shape))
#     print(original_data.head())
    original_data['orf'] = original_data['ORF'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    # Make sure everything translated ok
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data = original_data.loc[t,:]
    original_data['data'] = -1
    original_data.set_index('orf', inplace=True)
    original_data = original_data[['data']].copy()
    original_data = original_data.groupby(original_data.index).mean()
    
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[14]:


original_data = pd.concat(original_data_list, axis=1)


# In[16]:


original_data.index.name='orf'


# In[17]:


original_data[original_data.isnull()] = 0


# In[18]:


original_data.head()


# In[19]:


original_data.shape


# # Prepare the final dataset

# In[20]:


data = original_data.copy()


# In[21]:


dataset_ids = [497, 5004, 5005]
datasets = datasets.reindex(index=dataset_ids)


# In[22]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[23]:


data.head()


# ## Subset to the genes currently in SGD

# In[24]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[25]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[26]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[27]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[28]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

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



