#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 24035500
paper_name = 'shimada_gasser_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data1 = pd.read_excel('raw_data/Shimada_etal_BHS345_HIPHOP_Data.xlsx', sheet_name='HIP')
original_data2 = pd.read_excel('raw_data/Shimada_etal_BHS345_HIPHOP_Data.xlsx', sheet_name='HOP')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[7]:


orf_col = 'SYSTEMATIC_NAME'


# In[8]:


original_data1[orf_col] = original_data1[orf_col].astype(str)
original_data2[orf_col] = original_data2[orf_col].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data1[orf_col] = clean_orf(original_data1[orf_col])
original_data2[orf_col] = clean_orf(original_data2[orf_col])


# In[10]:


# Translate to ORFs 
original_data1[orf_col] = translate_sc(original_data1[orf_col], to='orf')
original_data2[orf_col] = translate_sc(original_data2[orf_col], to='orf')


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data1[orf_col])
print(original_data1.loc[~t,])


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data2[orf_col])
print(original_data2.loc[~t,])


# In[13]:


original_data1['data'] = original_data1['Z_SCORE']
original_data2['data'] = original_data2['Z_SCORE']


# In[14]:


original_data1.set_index(orf_col, inplace=True)
original_data2.set_index(orf_col, inplace=True)
original_data1.index.name='orf'
original_data2.index.name='orf'


# In[15]:


original_data1 = original_data1[['data']]
original_data2 = original_data2[['data']]


# In[16]:


original_data1 = original_data1.groupby(original_data1.index).mean()
original_data2 = original_data2.groupby(original_data2.index).mean()


# In[17]:


original_data = original_data1[['data']].join(original_data2[['data']], how='outer', lsuffix='_hip', rsuffix='_hop')


# In[18]:


original_data.head()


# # Prepare the final dataset

# In[19]:


data = original_data.copy()


# In[20]:


dataset_ids = [16608, 16609]
datasets = datasets.reindex(index=dataset_ids)


# In[21]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[22]:


data.head()


# ## Subset to the genes currently in SGD

# In[23]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[24]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[26]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[27]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[28]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[29]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[30]:


from IO.save_data_to_db3 import *


# In[31]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




