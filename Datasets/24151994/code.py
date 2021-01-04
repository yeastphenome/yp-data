#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 24151994
paper_name = 'marek_korona_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# ### Maximum growth rate (MGR)

# In[5]:


original_data1 = pd.read_excel('raw_data/evo12196-sup-0001-tables1.xls', sheet_name='Table S5', skiprows=1)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))


# In[7]:


original_data1['ORF'] = original_data1['ORF'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data1['ORF'] = clean_orf(original_data1['ORF'])


# In[9]:


# Translate to ORFs 
original_data1['ORF'] = translate_sc(original_data1['ORF'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['ORF'])
print(original_data1.loc[~t,])


# In[11]:


original_data1['data'] = original_data1[['Replica 1.1','Replica 2.1','Replica 3.1']].mean(axis=1)


# In[12]:


original_data1.set_index('ORF', inplace=True)
original_data1.index.name='orf'


# In[13]:


original_data1 = original_data1.groupby(original_data1.index).mean()


# In[14]:


original_data1.shape


# ### Maximum chronological lifespan (MLS)

# In[15]:


original_data2 = pd.read_excel('raw_data/evo12196-sup-0003-tables3.xls', sheet_name='3678 del', skiprows=1)


# In[16]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[17]:


original_data2['ORF'] = original_data2['ORF'].astype(str)


# In[18]:


# Eliminate all white spaces & capitalize
original_data2['ORF'] = clean_orf(original_data2['ORF'])


# In[19]:


# Translate to ORFs 
original_data2['ORF'] = translate_sc(original_data2['ORF'], to='orf')


# In[20]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['ORF'])
print(original_data2.loc[~t,])


# In[21]:


original_data2['data'] = original_data2[['Replica 1.1','Replica 2.1','Replica 3.1']].mean(axis=1)


# In[22]:


original_data2.set_index('ORF', inplace=True)
original_data2.index.name = 'orf'


# In[23]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[24]:


original_data2.shape


# ### Combine

# In[25]:


original_data = original_data1[['data']].join(original_data2[['data']], how='outer', lsuffix='_mgr', rsuffix='_mls')


# In[26]:


original_data.head()


# # Prepare the final dataset

# In[27]:


data = original_data.copy()


# In[28]:


dataset_ids = [16564, 16565]
datasets = datasets.reindex(index=dataset_ids)


# In[29]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[30]:


data.head()


# ## Subset to the genes currently in SGD

# In[31]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[32]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[33]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[34]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[35]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[36]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[37]:


from IO.save_data_to_db3 import *


# In[38]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




