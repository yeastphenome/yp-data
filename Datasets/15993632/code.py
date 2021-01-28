#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 15993632
paper_name = 'corbacho_hernandez_2005' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[11]:


original_data = pd.read_csv('raw_data/Corbacho_tables.txt', sep='\t', encoding = "ISO-8859-1", skiprows=2)


# In[12]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[13]:


original_data.head()


# In[16]:


original_data['orf'] = original_data['ORF '].astype(str)


# In[17]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[18]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[19]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[20]:


original_data = original_data.loc[t,:]


# In[23]:


original_data['data'] = original_data['Ldb group ']


# In[25]:


original_data.loc[original_data['data']=='1 (alpha)','data'] = '1'


# In[26]:


original_data['data'] = pd.to_numeric(original_data['data'], errors='coerce')


# In[28]:


original_data.set_index('orf', inplace=True)


# In[29]:


original_data = original_data[['data']].copy()


# In[30]:


original_data = original_data.groupby(original_data.index).mean()


# In[32]:


original_data['data'] = original_data['data']-4


# In[33]:


original_data.shape


# # Prepare the final dataset

# In[34]:


data = original_data.copy()


# In[35]:


dataset_ids = [582]
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




