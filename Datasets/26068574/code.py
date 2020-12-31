#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 26068574
paper_name = 'li_breeden_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/table_1.xlsx', sheet_name='Sheet1')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data['orfs'] = original_data['ORF ID'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[9]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[11]:


original_data = original_data.loc[t,]


# In[12]:


original_data['data'] = -1


# In[13]:


original_data.set_index('orfs', inplace=True)
original_data.index.name = 'orf'


# In[14]:


original_data = original_data[['data']].copy()


# In[15]:


original_data = original_data.groupby(original_data.index).mean()


# In[16]:


original_data.shape


# # Load & process tested strains

# In[17]:


tested = pd.read_excel('raw_data/library screen.xlsx', sheet_name='Sheet1')


# In[18]:


tested['orfs'] = tested['ORF name'].astype(str)
tested['orfs'] = clean_orf(tested['orfs'])


# In[19]:


tested.loc[tested['orfs']=='YLR287-A','orfs'] = 'YLR287C-A'


# In[20]:


tested['orfs'] = translate_sc(tested['orfs'], to='orf')


# In[21]:


# Make sure everything translated ok
t = looks_like_orf(tested['orfs'])
print(tested.loc[~t,])


# In[22]:


tested = tested.loc[t,]


# In[23]:


tested_orfs = np.unique(tested['orfs'].values)


# In[24]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]


# In[25]:


print(missing)


# In[26]:


# Adding the 2 missing hits to the list of tested strains
tested_orfs = np.array(list(tested_orfs) + missing)


# In[27]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[30]:


data = original_data.copy()


# In[31]:


dataset_ids = [16601]
datasets = datasets.reindex(index=dataset_ids)


# In[32]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[33]:


data.head()


# ## Subset to the genes currently in SGD

# In[34]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[35]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[36]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[37]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[38]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

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




