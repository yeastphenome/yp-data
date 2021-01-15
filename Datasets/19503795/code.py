#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 19503795
paper_name = 'westmoreland_bennett_2009' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data1 = pd.read_excel('raw_data/journal.pone.0005830.s001.XLS', sheet_name='Sheet1')
original_data2 = pd.read_excel('raw_data/journal.pone.0005830.s002.XLS', sheet_name='Sheet1')


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
t1 = looks_like_orf(original_data1['ORF'])
print(original_data1.loc[~t1,])


# In[11]:


# Make sure everything translated ok
t2 = looks_like_orf(original_data2['ORF'])
print(original_data2.loc[~t2,])


# In[12]:


original_data1 = original_data1.loc[t1,:]
original_data2 = original_data2.loc[t2,:]


# In[13]:


original_data1.set_index('ORF', inplace=True)
original_data2.set_index('ORF', inplace=True)

original_data1.index.name='orf'
original_data2.index.name='orf'


# In[14]:


original_data1['data'] = original_data1['DoxS (2n)'].astype(float)
original_data2['data'] = original_data2['DoxS (2n)'].astype(float)


# In[15]:


# "3" denotes the most sensitive strains; in addition, shifting all values by 2 to make them more extreme and unify with data2 (less sensitive phenotypes)
original_data1['data'] = -original_data1['data']-2
original_data2['data'] = -original_data2['data']


# In[16]:


original_data1 = original_data1[['data']].copy()
original_data2 = original_data2[['data']].copy()


# In[17]:


original_data1 = original_data1.groupby(original_data1.index).mean()
original_data2 = original_data2.groupby(original_data2.index).mean()


# In[18]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[19]:


original_data.head()


# In[21]:


original_data['data'] = original_data.mean(axis=1)
original_data = original_data[['data']].copy()


# In[22]:


original_data.shape


# # Prepare the final dataset

# In[23]:


data = original_data.copy()


# In[24]:


dataset_ids = [16454]
datasets = datasets.reindex(index=dataset_ids)


# In[25]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[26]:


data.head()


# ## Subset to the genes currently in SGD

# In[27]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[28]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[29]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[30]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[31]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[32]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[33]:


from IO.save_data_to_db3 import *


# In[34]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




