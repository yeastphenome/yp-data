#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 32054832
paper_name = 'chen_petranovic_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data1 = pd.read_excel('raw_data/complete score data for SGA.xlsx', sheet_name='nonessential repeat 1')
original_data2 = pd.read_excel('raw_data/complete score data for SGA.xlsx', sheet_name='nonessential repeat 2')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[7]:


original_data1['orfs'] = original_data1['Array ORF'].astype(str)
original_data2['orfs'] = original_data2['Array ORF'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data1['orfs'] = clean_orf(original_data1['orfs'])
original_data2['orfs'] = clean_orf(original_data2['orfs'])


# In[9]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['orfs'], to='orf')
original_data2['orfs'] = translate_sc(original_data2['orfs'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[12]:


data1 = original_data1[['orfs','Score']].groupby('orfs').mean()
data2 = original_data2[['orfs','Score']].groupby('orfs').mean()


# In[13]:


data = data1.join(data2, lsuffix='_1', rsuffix='_2')


# In[14]:


data['data'] = data[['Score_1','Score_2']].mean(axis=1)


# In[15]:


data = data.drop(columns=['Score_1','Score_2'])


# In[17]:


data.index.name='orf'


# # Prepare the final dataset

# In[18]:


dataset_ids = [16408]
datasets = datasets.reindex(index=dataset_ids)


# In[19]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[20]:


print('Final data dimensions: %d x %d' % (data.shape))


# ## Subset to the genes currently in SGD

# In[21]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[22]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# # Normalize

# In[23]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[24]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[25]:


data_norm[data.isnull()] = np.nan


# In[26]:


data_all = data.join(data_norm)


# In[27]:


data_all.head()


# # Print out

# In[28]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[29]:


from IO.save_data_to_db3 import *


# In[30]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




