#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 19638689
paper_name = 'auesukaree_harashima_2009' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


files = ['ethanol_sensitivity_hits.txt','methanol_sensitivity_hits.txt',
         'propanol_sensitivity_hits.txt','nacl_sensitivity_hits.txt',
         'h2o2_sensitivity_hits.txt','heat_sensitivity_hits.txt']


# In[10]:


original_data_list = []
for f in files:
    original_data = pd.read_csv('raw_data/'+f, header=None, sep='\t')
    print('Original data dimensions: %d x %d' % (original_data.shape))
#     print(original_data.head())
    original_data['gene'] = original_data[0].astype(str)
    original_data['gene'] = clean_genename(original_data['gene'])
    original_data['orf'] = translate_sc(original_data['gene'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data['data'] = -1
    original_data.set_index('orf', inplace=True)
    original_data = original_data[['data']].copy()
    original_data = original_data.groupby(original_data.index).mean()
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[11]:


original_data = pd.concat(original_data_list, axis=1)


# In[12]:


original_data.head()


# In[13]:


original_data.shape


# In[14]:


original_data[original_data.isnull()] = 0


# In[32]:


original_data.index.name = 'orf'


# # Load & process tested strains

# In[17]:


tested = pd.read_excel('raw_data/Mat alpha_KOset list.xlsx', sheet_name='Sheet1', skiprows=1)


# In[18]:


tested.head()


# In[19]:


tested['orf'] = tested['ORF name'].astype(str)


# In[20]:


tested['orf'] = clean_orf(tested['orf'])


# In[21]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[22]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[23]:


tested = tested.loc[t,:]


# In[24]:


tested_orfs = tested['orf'].unique()


# In[25]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[26]:


tested_orfs = list(tested_orfs) + missing


# In[27]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[33]:


data = original_data.copy()


# In[34]:


dataset_ids = [162, 432, 433, 434, 435, 436]
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


data_norm = normalize_phenotypic_scores(data, has_tested=True)


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




