#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 28592509
paper_name = 'acton_giaever_2017' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/rsob160330_si_003.xlsx', sheet_name='Additional File 3')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[12]:


# Drop the "_" annotations
original_data['ORF'] = original_data['ORF'].apply(lambda x: x.split('_')[0])


# In[13]:


original_data['orf'] = original_data['ORF'].astype(str)


# In[14]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[15]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[16]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[17]:


data_cols = [col for col in original_data.columns.values if col.startswith('YKO')]
data_cols


# In[18]:


original_data.set_index('orf', inplace=True)


# In[19]:


original_data = original_data[data_cols].copy()


# In[20]:


original_data = original_data.groupby(original_data.index).mean()


# In[21]:


original_data.shape


# In[35]:


# Switching sign to follow convention (sensitive strains = negative values)
original_data = -original_data


# In[36]:


original_data.head()


# # Prepare the final dataset

# In[37]:


data = original_data.copy()


# In[38]:


dataset_ids = [22009, 22010, 22011, 22012, 22013, 22014, 22015, 22016, 22017, 22018]
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


data_norm = normalize_phenotypic_scores(data, has_tested=True)


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




