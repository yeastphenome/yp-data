#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 22959270
paper_name = 'zhang_lobachev_2012' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, names=['genes'])


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


# In[11]:


original_data['data'] = 1


# In[12]:


original_data.set_index('orfs', inplace=True)


# In[13]:


original_data.index.name = 'orf'


# In[14]:


original_data = original_data.groupby(original_data.index).mean()


# In[15]:


original_data.shape


# # Load & process tested strains

# In[16]:


tested = pd.read_excel('raw_data/mat_a_101501.xlsx', sheet_name='mat_a_101501', skiprows=1)


# In[17]:


tested['orf'] = tested['ORF name'].astype(str)


# In[18]:


tested['orf'] = clean_orf(tested['orf'])


# In[19]:


tested.loc[tested['orf'] == 'YLR287-A','orf'] = 'YLR287C-A'


# In[20]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[21]:


t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[22]:


tested = tested.loc[t,:]


# In[23]:


tested_orfs = tested['orf'].unique()


# In[24]:


# Test if all hits are present in tested
missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
print(missing)


# In[25]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[26]:


data = original_data.copy()


# In[27]:


dataset_ids = [16402]


# In[28]:


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




