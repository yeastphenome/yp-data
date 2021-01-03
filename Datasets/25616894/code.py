#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 25616894
paper_name = 'junne_hoepfner_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


sheet_names = ['Decatransin Cpd 1 Exp 1','Decatransin Cpd 2 Exp 2','Cotransin Cpd 2','Cotansin Cpd 3 (HUN-7293)']


# In[7]:


original_data1_list = []
original_data2_list = []
for s in sheet_names:
    original_data = pd.read_excel('raw_data/Junnet_et_al_HIPHOPrawdata.xlsx', sheet_name=s)
    print('Original data dimensions: %d x %d' % (original_data.shape))
    original_data['orf'] = original_data['SYSTEMATIC_NAME'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    
    # Make sure everything translated ok
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    
    original_data['data'] = original_data['Z_SCORE']
    original_data.set_index('orf', inplace=True)
    
    original_data1 = original_data.loc[original_data['EXPERIMENT_TYPE']=='HIP',:].copy()
    original_data2 = original_data.loc[original_data['EXPERIMENT_TYPE']=='HOP',:].copy()
    
    original_data1 = original_data1[['data']].copy()
    original_data2 = original_data2[['data']].copy()
    
    original_data1 = original_data1.groupby(original_data1.index).mean()
    original_data2 = original_data2.groupby(original_data2.index).mean()
    
    print(original_data1.shape)
    print(original_data2.shape)
    
    original_data1_list.append(original_data1)
    original_data2_list.append(original_data2)


# In[17]:


original_data1 = pd.concat(original_data1_list, axis=1)
original_data2 = pd.concat(original_data2_list, axis=1)


# In[18]:


original_data1.columns = np.arange(4)
original_data1[4] = original_data1[[0,1]].mean(axis=1)
original_data1.drop(columns=[0,1], inplace=True)
original_data1.head()


# In[20]:


original_data2.columns = np.arange(4)
original_data2[4] = original_data2[[0,1]].mean(axis=1)
original_data2.drop(columns=[0,1], inplace=True)
original_data2.head()


# In[21]:


original_data = original_data1.join(original_data2, lsuffix='_1', rsuffix='_2')


# In[24]:


original_data.index.name='orf'


# In[25]:


original_data.head()


# # Prepare the final dataset

# In[27]:


data = original_data.copy()


# In[28]:


dataset_ids = [5265,5267,773,5264,5266,772]
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




