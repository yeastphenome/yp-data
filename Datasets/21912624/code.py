#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 21912624
paper_name = 'north_vulpe_2011' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


files = ['s007.xlsx','s008.xlsx','s009.xlsx']


# In[11]:


original_data_list = []
for f in files:
    original_data = pd.read_excel('raw_data/' + f, sheet_name='Sheet1', skiprows=1)
    print('Original data dimensions: %d x %d' % (original_data.shape))
    original_data['orf'] = original_data['ORF'].astype(str)
    # Eliminate all white spaces & capitalize
    original_data['orf'] = clean_orf(original_data['orf'])
    # Translate to ORFs 
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    # Make sure everything translated ok
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data = original_data.loc[t,:]
    
    original_data.set_index('orf', inplace=True)
    original_data = original_data.iloc[:,2:8].copy()
    
    for c in original_data.columns:
        original_data[c] = pd.to_numeric(original_data[c], errors='coerce')
        
    original_data = original_data.groupby(original_data.index).mean()
    original_data_list.append(original_data)


# In[19]:


original_data = pd.concat(original_data_list, axis=1)


# In[20]:


original_data.head()


# In[21]:


original_data.shape


# In[22]:


original_data.index.name='orf'


# In[27]:


original_data[original_data.isnull()] = 0


# # Prepare the final dataset

# In[28]:


data = original_data.copy()


# In[29]:


dataset_ids = [513, 1330, 1331, 775, 1332, 1333,
               1334, 1335, 1336, 1337, 1338, 1339, 
               1340, 1342, 1341, 1343, 1345, 1344]
datasets = datasets.reindex(index=dataset_ids)


# In[30]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[31]:


data.head()


# ## Subset to the genes currently in SGD

# In[32]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[33]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[35]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


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




