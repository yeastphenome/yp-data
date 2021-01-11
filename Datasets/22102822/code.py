#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 22102822
paper_name = 'berry_gasch_2011' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/pgen.1002353.s009.xlsx', sheet_name='Hom-Het COMPILATION')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data['orf'] = original_data['YORF'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[10]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[12]:


original_data.set_index('orf', inplace=True)


# In[13]:


original_data.drop(columns=['YORF'], inplace=True)


# In[14]:


for c in original_data.columns:
    original_data[c] = pd.to_numeric(original_data[c], errors='coerce')


# ## Data processing

# In[15]:


import re


# In[16]:


regex_list = ['Sample1 vs Sample0 Array','Sample2 vs Sample1 Array',
              'DTT[0-9] Sample3 vs Sample1 Array','NaCl[0-9] Sample3 vs Sample1 Array','HS[0-9] Sample3 vs Sample1 Array',
              'DTT[0-9] Sample4[A-Z]? vs Sample3 Array','NaCl[0-9] Sample4[A-Z]? vs Sample3 Array','HS[0-9] Sample4[A-Z]? vs Sample3 Array',
              'Sample1 vs Sample0 [(DN)(UP)]','Sample2 vs Sample1 [(DN)(UP)]','DTT[0-9] Sample3 vs Sample1 [(DN)(UP)]',
              'NaCl[0-9] Sample3 vs Sample1 [(DN)(UP)]','HS[0-9] Sample3 vs Sample1 [(DN)(UP)]','TM[0-9] Sample3 vs Sample1 [(DN)(UP)]',
              'DTT[0-9] Sample4[A-Z]? vs Sample3 [(DN)(UP)]','NaCl[0-9] Sample4[A-Z]? vs Sample3 [(DN)(UP)]',
              'HS[0-9] Sample4[A-Z]? vs Sample3 [(DN)(UP)]','TM[0-9] Sample4[A-Z]? vs Sample3 [(DN)(UP)]']


# In[17]:


original_data_list = []
for r in regex_list:
    data_cols = [c for c in original_data.columns if bool(re.search(r, c))]
    t = original_data[data_cols].mean(axis=1)
    original_data_list.append(t)


# In[18]:


original_data = pd.concat(original_data_list, axis=1)


# In[19]:


original_data = original_data.groupby(original_data.index).mean()


# In[20]:


### Only keep Hom strains using current Hom collection from Open Biosystems
hom = pd.read_excel('extras/Homozygous_diploid_obs_v7.0.xlsx', sheet_name='DATA')
hom.head()


# In[21]:


hom['orf'] = hom['ORF'].astype(str)


# In[22]:


hom['orf'] = clean_orf(hom['orf'])


# In[23]:


hom['orf'] = translate_sc(hom['orf'], to='orf')


# In[24]:


# Make sure everything translated ok
t = looks_like_orf(hom['orf'])
print(hom.loc[~t,])


# In[25]:


hom_orfs = hom['orf'].unique()


# In[26]:


hom_orfs_in_data = [orf for orf in hom_orfs if orf in original_data.index.values]
len(hom_orfs_in_data)


# In[27]:


original_data = original_data.reindex(index=hom_orfs_in_data, fill_value=np.nan)


# # Prepare the final dataset

# In[28]:


data = original_data.copy()


# In[29]:


dataset_ids = [758, 759, 761, 760, 762, 590, 588, 589, 5395, 5396, 5397, 5398, 5399, 5400, 5401, 5402, 5403, 5404]
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

# In[34]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[35]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[36]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[37]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[38]:


from IO.save_data_to_db3 import *


# In[39]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




