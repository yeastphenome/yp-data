#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 15388447
paper_name = 'markovich_osherov_2004' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[11]:


original_data = pd.read_csv('raw_data/hits.txt', sep='\t', header=None)


# In[12]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[13]:


original_data.head()


# In[14]:


sensitive = original_data.loc[0,0].split(',')
resistant = original_data.loc[1,0].split(',')


# In[15]:


original_data1 = pd.DataFrame(data={'gene': sensitive,'data': -1})
original_data2 = pd.DataFrame(data={'gene': resistant,'data': 1})


# In[22]:


original_data = pd.concat([original_data1, original_data2], axis=0)


# In[23]:


original_data['gene'] = original_data['gene'].astype(str)


# In[24]:


original_data['gene'] = original_data['gene'].apply(lambda x: x.split('/')[0])


# In[25]:


# Eliminate all white spaces & capitalize
original_data['gene'] = clean_genename(original_data['gene'])


# In[26]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['gene'].values, to='orf')


# In[27]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[28]:


original_data = original_data.loc[t,:]


# In[29]:


original_data.set_index('orf', inplace=True)


# In[30]:


original_data = original_data[['data']].copy()


# In[31]:


original_data = original_data.groupby(original_data.index).mean()


# In[32]:


original_data.shape


# # Prepare the final dataset

# In[33]:


data = original_data.copy()


# In[34]:


dataset_ids = [189]
datasets = datasets.reindex(index=dataset_ids)


# In[35]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[36]:


data.head()


# ## Subset to the genes currently in SGD

# In[37]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[38]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[39]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[40]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[41]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[42]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[43]:


from IO.save_data_to_db3 import *


# In[44]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:



