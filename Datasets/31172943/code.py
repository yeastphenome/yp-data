#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 31172943
paper_name = 'dederer_lamberg_2019' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data1 = pd.read_excel('raw_data/elife-45506-supp1-v2.xlsx', sheet_name='tFT-Pex15delta30')
original_data2 = pd.read_excel('raw_data/elife-45506-supp1-v2.xlsx', sheet_name='tFT-Tom5')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[7]:


original_data1['ORF'] = original_data1['ORF'].astype(str)
original_data2['ORF'] = original_data2['ORF'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data1['ORF'] = clean_orf(original_data1['ORF'])
original_data2['ORF'] = clean_orf(original_data2['ORF'])


# In[9]:


# Translate to ORFs 
original_data1['ORF'] = translate_sc(original_data1['ORF'], to='orf')
original_data2['ORF'] = translate_sc(original_data2['ORF'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['ORF'])
print(original_data1.loc[~t,])


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['ORF'])
print(original_data2.loc[~t,])


# In[12]:


original_data1 = original_data1.groupby('ORF').mean()
original_data1.index.name='orf'


# In[13]:


original_data2 = original_data2.groupby('ORF').mean()
original_data2.index.name='orf'


# In[14]:


data = original_data1[['ratio.mean']].join(original_data2[['ratio.mean']], how='outer', lsuffix='_1', rsuffix='_2')


# # Prepare the final dataset

# In[15]:


dataset_ids = [16549,16550]
datasets = datasets.reindex(index=dataset_ids)


# In[16]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[17]:


data.head()


# ## Subset to the genes currently in SGD

# In[18]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[19]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[20]:


data.head()


# # Normalize

# In[21]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[22]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[23]:


data_norm[data.isnull()] = np.nan


# In[24]:


data_all = data.join(data_norm)


# In[25]:


data_all.head()


# # Print out

# In[26]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[27]:


from IO.save_data_to_db3 import *


# In[28]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




