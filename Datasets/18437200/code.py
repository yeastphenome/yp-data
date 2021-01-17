#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 18437200
paper_name = 'jin_freedman_2008' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


sheets = ['EC10','EC50']


# In[10]:


data_cols


# In[12]:


original_data_list = []
for s in sheets:
    original_data = pd.read_excel('raw_data/journal.pgen.1000053.s005.xlsx', sheet_name=s, skiprows=3)
    print('Original data dimensions: %d x %d' % (original_data.shape))
#     print(original_data.head())
    if 'NAME' in original_data.columns:
        original_data['gene'] = original_data['NAME'].astype(str)
    else:
        original_data['gene'] = original_data['Gene'].astype(str)
        
    original_data['gene'] = clean_genename(original_data['gene'])
    original_data['orf'] = translate_sc(original_data['gene'], to='orf')
    # Make sure everything translated ok
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data.set_index('orf', inplace=True)
    
    data_cols = original_data.columns.values[2:9]
    original_data = original_data[data_cols].astype(float)
    
    original_data = original_data.groupby(original_data.index).mean()
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[13]:


original_data = pd.concat(original_data_list, axis=1)


# In[16]:


original_data.index.name='orf'


# In[17]:


original_data = (1 / original_data) - 1


# In[30]:


original_data[original_data.isnull()] = 0


# # Prepare the final dataset

# In[31]:


data = original_data.copy()


# In[32]:


dataset_ids = [11772, 11773, 11771, 11774, 11775, 11776, 11777, 1311, 1312, 1313, 1314, 1310, 1315, 1316]
datasets = datasets.reindex(index=dataset_ids)


# In[33]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[34]:


data.head()


# ## Subset to the genes currently in SGD

# In[35]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[36]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[37]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[38]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[39]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[40]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[41]:


from IO.save_data_to_db3 import *


# In[42]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




