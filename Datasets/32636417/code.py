#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 32636417
paper_name = 'fu_hoepfner_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


sheets = ['Raw Data HIP HOP Exp.1 Fig 2ab','Raw Data HIP HOP Exp.2 Fig 2ab']


# In[6]:


original_data1_list = []
original_data2_list = []

for s in sheets:
    original_data = pd.read_excel('raw_data/SourceData_TableFormats.xlsx', sheet_name=s)
    print('Original data dimensions: %d x %d' % (original_data.shape))
    original_data.head()
    original_data['orf'] = original_data['SYSTEMATIC_NAME'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data['data'] = pd.to_numeric(original_data['SCORE'], errors='coerce')
    original_data.set_index('orf', inplace=True)
    
    original_data1 = original_data.loc[original_data['EXPERIMENT_TYPE']=='HOP'].copy()
    original_data2 = original_data.loc[original_data['EXPERIMENT_TYPE']=='HIP'].copy()
    
    original_data1 = original_data1[['data']].copy()
    original_data2 = original_data2[['data']].copy()
    
    original_data1 = original_data1.groupby(original_data1.index).mean()
    original_data2 = original_data2.groupby(original_data2.index).mean()
    
    print(original_data1.shape)
    print(original_data2.shape)
    
    original_data1_list.append(original_data1)
    original_data2_list.append(original_data2)


# In[7]:


original_data1 = pd.concat(original_data1_list, axis=0)
original_data1 = original_data1.groupby(original_data1.index).mean()
original_data1.shape


# In[8]:


original_data2 = pd.concat(original_data2_list, axis=0)
original_data2 = original_data2.groupby(original_data2.index).mean()
original_data2.shape


# In[9]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[10]:


original_data.shape


# # Prepare the final dataset

# In[11]:


data = original_data.copy()


# In[12]:


dataset_ids = [21880, 21879]
datasets = datasets.reindex(index=dataset_ids)


# In[13]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[14]:


data.head()


# ## Subset to the genes currently in SGD

# In[15]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[16]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[17]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[18]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[19]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[20]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[21]:


from IO.save_data_to_db3 import *


# In[22]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




