#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 22704625
paper_name = 'hoepfner_winzeler_2012' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[6]:


sheets = ['Cladosporin Exp 1','Cladosporin Exp 2','Cladosporin Exp 3']


# In[9]:


original_data_list1 = []
original_data_list2 = []
for s in sheets:
    original_data = pd.read_excel('raw_data/Hoepfner_et_al_HIPHOPrawdate.xlsx', sheet_name=s)
    print('Original data dimensions: %d x %d' % (original_data.shape))
    
    original_data['orf'] = original_data['SYSTEMATIC_NAME'].astype(str)
    # Eliminate all white spaces & capitalize
    original_data['orf'] = clean_orf(original_data['orf'])
    # Translate to ORFs 
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    # Make sure everything translated ok
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    
    original_data1 = original_data.loc[original_data['EXPERIMENT_TYPE']=='HOP',:].copy()
    original_data2 = original_data.loc[original_data['EXPERIMENT_TYPE']=='HIP',:].copy()
    
    original_data1['data'] = -original_data1['SENSITIVITY SCORE']
    original_data2['data'] = -original_data2['SENSITIVITY SCORE']    
    
    original_data1.set_index('orf', inplace=True)
    original_data2.set_index('orf', inplace=True)
    
    original_data1 = original_data1[['data']].copy()
    original_data2 = original_data2[['data']].copy()
    
    original_data1 = original_data1.groupby(original_data1.index).mean()
    original_data2 = original_data2.groupby(original_data2.index).mean()
    
    original_data_list1.append(original_data1)
    original_data_list2.append(original_data2)


# In[20]:


original_data1 = pd.concat(original_data_list1, axis=1)
original_data2 = pd.concat(original_data_list2, axis=1)


# In[21]:


original_data1['avg'] = original_data1.mean(axis=1)
original_data2['avg'] = original_data2.mean(axis=1)


# In[22]:


original_data1 = original_data1[['avg']].copy()
original_data1.columns = ['data']


# In[23]:


original_data2 = original_data2[['avg']].copy()
original_data2.columns = ['data']


# In[24]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[27]:


original_data.index.name = 'orf'


# In[28]:


original_data.shape


# # Prepare the final dataset

# In[29]:


data = original_data.copy()


# In[30]:


dataset_ids = [5251, 5252]
datasets = datasets.reindex(index=dataset_ids)


# In[31]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


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

data.head()


# # Normalize

# In[35]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[36]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[37]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[38]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[39]:


from IO.save_data_to_db3 import *


# In[40]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




