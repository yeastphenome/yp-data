#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 14676322
paper_name = 'warringer_blomberg_2003' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[8]:


original_data = pd.read_excel('raw_data/LPI NaCl.xlsx', sheet_name='LPI', skiprows=3)


# In[9]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[10]:


original_data.head()


# In[11]:


original_data['orf'] = original_data['Unnamed: 0'].astype(str)


# In[12]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[13]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[14]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[15]:


original_data = original_data.loc[t,:]


# In[16]:


original_data.set_index('orf', inplace=True)


# In[17]:


original_data = original_data[['Adaptation','GrowthRate','Stationary Phase']].apply(pd.to_numeric, axis=1, errors='coerce')


# In[18]:


original_data = original_data.groupby(original_data.index).mean()


# In[19]:


original_data.shape


# # Load control data

# In[20]:


original_data2 = pd.read_excel('raw_data/LSC Reference.xlsx', sheet_name='LSC', skiprows=3)


# In[21]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[22]:


original_data2.head()


# In[23]:


original_data2['orf'] = original_data2['Unnamed: 0'].astype(str)


# In[24]:


# Eliminate all white spaces & capitalize
original_data2['orf'] = clean_orf(original_data2['orf'])


# In[25]:


# Translate to ORFs 
original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')


# In[26]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[27]:


original_data2 = original_data2.loc[t,:]


# In[28]:


original_data2.set_index('orf', inplace=True)


# In[31]:


original_data2 = original_data2[['Average','Average.1','Average.2']].apply(pd.to_numeric, axis=1, errors='coerce')


# In[32]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[33]:


original_data2.shape


# In[34]:


original_data = original_data.join(original_data2, how='outer')


# # Prepare the final dataset

# In[35]:


data = original_data.copy()


# In[36]:


dataset_ids = [50, 49, 51, 16184, 16183, 16185]
datasets = datasets.reindex(index=dataset_ids)


# In[37]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[38]:


data.head()


# ## Subset to the genes currently in SGD

# In[39]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[40]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[41]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[42]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[43]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[44]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[45]:


from IO.save_data_to_db3 import *


# In[46]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




