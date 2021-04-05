#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 18192430
paper_name = 'linderholm_bisson_2008' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[12]:


original_data = pd.read_csv('raw_data/table3.txt', sep='\t')


# In[13]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[14]:


original_data.head()


# In[15]:


original_data['genes'] = original_data['Gene(s)b'].apply(lambda x: x.split(', '))
genes = [x for _list in original_data['genes'] for x in _list ]


# In[17]:


original_data = pd.DataFrame(data={'genes': genes, 'data': 1})


# In[21]:


original_data2 = pd.read_csv('raw_data/white_colonies.txt', sep='\t', header=None, names=['genes','data'])
original_data2.head()


# In[22]:


original_data = pd.concat([original_data, original_data2], axis=0)


# In[28]:


original_data = original_data.reset_index()


# In[29]:


original_data['genes'] = original_data['genes'].astype(str)


# In[30]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[31]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['genes'], to='orf')


# In[32]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[33]:


original_data.set_index('orf', inplace=True)


# In[34]:


original_data = original_data[['data']].copy()


# In[35]:


original_data = original_data.groupby(original_data.index).mean()


# In[36]:


original_data.shape


# # Prepare the final dataset

# In[37]:


data = original_data.copy()


# In[38]:


dataset_ids = [16706]
datasets = datasets.reindex(index=dataset_ids)


# In[39]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[40]:


data.head()


# ## Subset to the genes currently in SGD

# In[41]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[42]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[43]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[44]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[45]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[46]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[47]:


from IO.save_data_to_db3 import *


# In[48]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




