#!/usr/bin/env python
# coding: utf-8

# In[3]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[4]:


paper_pmid = 20520766
paper_name = 'emadi_vuica_ross_2010' 


# In[5]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[6]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[7]:


original_data = pd.read_excel('raw_data/pone.0010846.s005.xls', sheet_name='Sheet1', skiprows=1)


# In[8]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[9]:


original_data.head()


# In[20]:


original_data1 = original_data.iloc[0:5,:].copy()
original_data2 = original_data.iloc[5:8,:].copy()
original_data3 = original_data.iloc[8:13,:].copy()


# In[21]:


original_data_list = [original_data1, original_data2, original_data3]


# In[23]:


original_data_list2 = []
for original_data in original_data_list:
    original_data['orf'] = original_data['ORF'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data = original_data.loc[t,:]
    
    original_data['data'] = original_data['Score'].apply(lambda x: int(x[0]))
    original_data.set_index('orf', inplace=True)
    original_data = original_data[['data']].copy()
    original_data = original_data.groupby(original_data.index).mean()
    
    print(original_data.shape)
    
    original_data_list2.append(original_data)


# In[24]:


original_data = pd.concat(original_data_list2, axis=1)


# In[26]:


original_data[original_data.isnull()] = 0


# In[33]:


original_data.index.name = 'orf'


# # Prepare the final dataset

# In[34]:


data = original_data.copy()


# In[35]:


dataset_ids = [16632, 16633, 16631]
datasets = datasets.reindex(index=dataset_ids)


# In[36]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[37]:


data.head()


# ## Subset to the genes currently in SGD

# In[38]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[39]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[40]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[41]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[42]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[43]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[44]:


from IO.save_data_to_db3 import *


# In[45]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




