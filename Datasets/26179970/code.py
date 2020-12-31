#!/usr/bin/env python
# coding: utf-8

# In[12]:


get_ipython().run_line_magic('run', '../yp_utils.py')
from functools import reduce


# # Initial setup

# In[2]:


paper_pmid = 26179970
paper_name = 'krastel_hoepfner_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[15]:


sheet_names = ['Nannocystin Exp 1','Nannocystin Exp 2']


# In[16]:


original_data_list1 = []
original_data_list2 = []

for s in sheet_names:
    original_data = pd.read_excel('raw_data/Krastel_et_al_HIPHOPrawdata.xlsx', sheet_name=s)
    print('Original data dimensions: %d x %d' % (original_data.shape))
    original_data['orf'] = original_data['SYSTEMATIC_NAME'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    # Make sure everything translated ok
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data.set_index('orf', inplace=True)
    original_data['data'] = original_data['SENSITIVITY SCORE']
    
    original_data1 = original_data.loc[original_data['EXPERIMENT_TYPE']=='HOP',:].copy()
    original_data1 = original_data1[['data']].copy()

    original_data2 = original_data.loc[original_data['EXPERIMENT_TYPE']=='HIP',:].copy()
    original_data2 = original_data2[['data']].copy()
    
    original_data1 = original_data1.groupby(original_data1.index).mean()
    print(original_data1.shape)
    original_data_list1.append(original_data1)
    
    original_data2 = original_data2.groupby(original_data2.index).mean()
    print(original_data2.shape)
    
    original_data_list2.append(original_data2)


# In[17]:


original_data1 = reduce(lambda x, y: pd.merge(x, y, on='orf', how='outer'), original_data_list1)


# In[19]:


original_data2 = reduce(lambda x, y: pd.merge(x, y, on='orf', how='outer'), original_data_list2)


# In[20]:


original_data1['data'] = original_data1.mean(axis=1)
original_data2['data'] = original_data2.mean(axis=1)


# In[23]:


original_data = original_data1[['data']].join(original_data2[['data']], how='outer', lsuffix='_1', rsuffix='_2')


# In[24]:


original_data.head()


# # Prepare the final dataset

# In[25]:


data = original_data.copy()


# In[26]:


dataset_ids = [5253, 5254]
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


data_norm = normalize_phenotypic_scores(data, has_tested=True)


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




