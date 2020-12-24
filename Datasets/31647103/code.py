#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 31647103
paper_name = 'schmidt_hombauer_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data1 = pd.read_csv('raw_data/CAN1.txt', header=None, names=['genes'])
original_data2 = pd.read_csv('raw_data/lys2-10A.txt', header=None, names=['genes'])


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[7]:


original_data1['genes'] = original_data1['genes'].astype(str)
original_data2['genes'] = original_data2['genes'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data1['genes'] = clean_genename(original_data1['genes'])
original_data2['genes'] = clean_genename(original_data2['genes'])


# In[9]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['genes'], to='orf')
original_data2['orfs'] = translate_sc(original_data2['genes'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[12]:


original_data1['data'] = 1
original_data2['data'] = 1


# In[13]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)


# In[14]:


original_data1.index.name='orf'
original_data2.index.name='orf'


# # Load & process tested strains

# In[15]:


tested = pd.read_excel('raw_data/transomic collection.xlsx', sheet_name='list with names')


# In[16]:


tested['SystematicName'] = tested['SystematicName'].astype(str)


# In[17]:


tested['SystematicName'] = clean_orf(tested['SystematicName'])


# In[18]:


tested['SystematicName'] = translate_sc(tested['SystematicName'], to='orf')


# In[19]:


t = looks_like_orf(tested['SystematicName'])
print(tested.loc[~t,])


# In[20]:


tested_orfs = tested['SystematicName'].unique()


# In[21]:


missing = [orf for orf in original_data1.index.values if orf not in tested_orfs]
missing


# In[22]:


missing = [orf for orf in original_data2.index.values if orf not in tested_orfs]
missing


# # Prepare the final dataset

# In[23]:


dataset_ids = [16441, 16442]
datasets = datasets.reindex(index=dataset_ids)


# In[24]:


data = pd.DataFrame(index=tested_orfs, columns=dataset_ids, data=0)


# In[25]:


data.head()


# In[26]:


data.loc[original_data1.index, dataset_ids[0]] = original_data1['data'].values


# In[27]:


data.loc[original_data2.index, dataset_ids[1]] = original_data2['data'].values


# In[28]:


data = data.groupby(data.index.values).mean()


# In[31]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[30]:


print('Final data dimensions: %d x %d' % (data.shape))


# In[32]:


data.head()


# ## Subset to the genes currently in SGD

# In[33]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[34]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[35]:


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


# In[39]:


data_all = data.join(data_norm)


# In[40]:


data_all.head()


# # Print out

# In[41]:


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




