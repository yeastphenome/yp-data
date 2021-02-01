#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 12149442
paper_name = 'hanway_romesberg_2002' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


sheet_names = ['UV','MMS']


# In[9]:


original_data_list = []
for s in sheet_names:
    original_data = pd.read_excel('raw_data/uv_hits.xlsx', sheet_name=s)
    print('Original data dimensions: %d x %d' % (original_data.shape))
#     print(original_data.head())
    original_data['gene'] = original_data['ORF'].astype(str)
    original_data['gene'] = clean_orf(original_data['gene'])
    
    typo_fixes = {'NTG1D':'NTG1','TOS10': 'YGR153W'}
    original_data['gene'] = original_data['gene'].apply(lambda x: typo_fixes[x] if x in typo_fixes.keys() else x)
    
    original_data['orf'] = translate_sc(original_data['gene'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data = original_data.loc[t,:]
    
    original_data.set_index('orf', inplace=True)
    
    data_cols = [c for c in original_data.columns if c.startswith('RF')]
    original_data = original_data[data_cols].apply(pd.to_numeric, axis=1, errors='coerce')
    original_data = original_data.groupby(original_data.index).mean()
    
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[10]:


original_data = pd.concat(original_data_list, axis=1)


# In[12]:


original_data.index.name='orf'

# Converting to scale so that wt = 0 and negative values = sensitivity
original_data = original_data-1


# In[13]:


original_data.shape


# In[24]:


# Assuming all non-hits were tested in the other screen
original_data[original_data.isnull()] = 0


# # Prepare the final dataset

# In[25]:


data = original_data.copy()


# In[26]:


dataset_ids = [474, 11850, 11851, 11852, 11853, 11854, 11855, 11856, 4951, 4967]
datasets = datasets.reindex(index=dataset_ids)


# In[27]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[28]:


data.head()


# ## Subset to the genes currently in SGD

# In[29]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[30]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[31]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[32]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[33]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[34]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[35]:


from IO.save_data_to_db3 import *


# In[36]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




