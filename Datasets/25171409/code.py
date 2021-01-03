#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 25171409
paper_name = 'jarosz_lindquist_2014' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data1 = pd.read_excel('raw_data/mmc3.xls', sheet_name='Sheet1')
original_data2 = pd.read_excel('raw_data/mmc4.xlsx', sheet_name='Sheet1')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[7]:


original_data1['orf'] = original_data1['ORF name'].astype(str)
original_data2['orf'] = original_data2['ORF name'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data1['orf'] = clean_orf(original_data1['orf'])
original_data2['orf'] = clean_orf(original_data2['orf'])


# In[9]:


# Translate to ORFs 
original_data1['orf'] = translate_sc(original_data1['orf'], to='orf')
original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orf'])
print(original_data1.loc[~t,])


# In[11]:


original_data1.drop(index=original_data1.loc[~t].index, inplace=True)


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[13]:


original_data2.drop(index=original_data2.loc[~t].index, inplace=True)


# In[14]:


original_data1['data'] = -1
original_data2['data'] = 1


# In[15]:


original_data1.set_index('orf', inplace=True)
original_data2.set_index('orf', inplace=True)


# In[16]:


original_data1 = original_data1[['data']].copy()
original_data2 = original_data2[['data']].copy()


# In[17]:


original_data1 = original_data1.groupby(original_data1.index).mean()
original_data2 = original_data2.groupby(original_data2.index).mean()


# In[20]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[21]:


original_data.head()


# In[22]:


original_data['data'] = original_data.mean(axis=1)


# In[23]:


original_data.drop(columns=['data_1','data_2'], inplace=True)


# # Prepare the final dataset

# In[24]:


data = original_data.copy()


# In[25]:


dataset_ids = [16510]
datasets = datasets.reindex(index=dataset_ids)


# In[26]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[27]:


data.head()


# ## Subset to the genes currently in SGD

# In[28]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[29]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[30]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[31]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[32]:


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




