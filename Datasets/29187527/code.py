#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 29187527
paper_name = 'eisenberg_bord_bohnert_2018' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# ### Part 1

# In[5]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, names=['genes','data'], sep='\t')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data['genes'] = original_data['genes'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[9]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[13]:


original_data.set_index('orfs', inplace=True)
original_data.index.name='orf'


# ### Part 2

# In[14]:


original_data2 = pd.read_excel('raw_data/JCB_201704122_TableS1.xlsx', sheet_name='HITS from screens', skiprows=2, nrows=30)


# In[15]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[16]:


original_data2['ORF'] = original_data2['ORF'].astype(str)


# In[17]:


# Eliminate all white spaces & capitalize
original_data2['ORF'] = clean_orf(original_data2['ORF'])


# In[18]:


# Translate to ORFs 
original_data2['ORF'] = translate_sc(original_data2['ORF'], to='orf')


# In[19]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['ORF'])
print(original_data2.loc[~t,])


# In[20]:


original_data2['data'] = 1


# In[21]:


original_data2.set_index('ORF', inplace=True)
original_data2.index.name='orf'


# ### Merge

# In[22]:


data = original_data[['data']].join(original_data2[['data']], lsuffix='_lower', rsuffix='_higher', how='outer')


# In[23]:


data['data'] = data[['data_lower','data_higher']].mean(axis=1)


# In[24]:


data = data[['data']]


# ### Remove essential genes (if any) -- original data included deletions and DAmP strains

# In[25]:


ess = ~is_essential(data.index.values)


# In[26]:


data = data.loc[ess.values,:]


# In[27]:


data = data.groupby(data.index).mean()


# In[28]:


data.shape


# # Prepare the final dataset

# In[29]:


dataset_ids = [16438]
datasets = datasets.reindex(index=dataset_ids)


# In[30]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[31]:


data.head()


# ## Subset to the genes currently in SGD

# In[32]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[33]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[34]:


data.head()


# # Normalize

# In[35]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[36]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[37]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[38]:


data_all.head()


# # Print out

# In[39]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[40]:


from IO.save_data_to_db3 import *


# In[41]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




