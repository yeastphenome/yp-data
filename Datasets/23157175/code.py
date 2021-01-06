#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 23157175
paper_name = 'zhang_jiang_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data1 = pd.read_excel('raw_data/hits.xlsx', sheet_name='DMSO4')
original_data2 = pd.read_excel('raw_data/hits.xlsx', sheet_name='DMSO8')


# In[6]:


original_data1.columns = [x.strip() for x in original_data1.columns]
original_data2.columns = [x.strip() for x in original_data2.columns]


# In[7]:


original_data1.drop(index=original_data1.loc[original_data1['DMSO sensitivity'].isnull()].index, inplace=True)
original_data2.drop(index=original_data2.loc[original_data2['8% DMSO sensitivity'].isnull()].index, inplace=True)


# In[8]:


original_data1['DMSO sensitivity'] = original_data1['DMSO sensitivity'].apply(lambda x: len(x.strip()))
original_data2['8% DMSO sensitivity'] = original_data2['8% DMSO sensitivity'].apply(lambda x: len(x.strip()))


# In[9]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[10]:


original_data1['genes'] = original_data1['Gene'].astype(str)
original_data2['genes'] = original_data2['Gene'].astype(str)


# In[11]:


# Eliminate all white spaces & capitalize
original_data1['genes'] = clean_genename(original_data1['genes'])
original_data2['genes'] = clean_genename(original_data2['genes'])


# In[12]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['genes'].values, to='orf')
original_data2['orfs'] = translate_sc(original_data2['genes'].values, to='orf')


# In[13]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[14]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[15]:


original_data1['data'] = -original_data1['DMSO sensitivity']
original_data2['data'] = -original_data2['8% DMSO sensitivity']


# In[16]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)
original_data1.index.name='orf'
original_data2.index.name='orf'


# In[17]:


original_data1 = original_data1[['data']].copy()
original_data2 = original_data2[['data']].copy()


# In[18]:


original_data1 = original_data1.groupby(original_data1.index).mean()
original_data2 = original_data2.groupby(original_data2.index).mean()


# In[20]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[23]:


original_data[original_data.isnull()] = 0


# In[24]:


original_data.head()


# In[25]:


original_data.shape


# # Load & process tested strains

# In[26]:


tested = pd.read_excel('raw_data/DELETION LIBRARY.xlsx', sheet_name='DELETION LIBRARY', skiprows=1)


# In[27]:


tested.head()


# In[28]:


tested['orf'] = tested['ORF name'].astype(str)


# In[29]:


tested['orf'] = clean_orf(tested['orf'])


# In[30]:


tested.loc[tested['orf'] == 'YELOO1C','orf'] = 'YEL001C'


# In[31]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[32]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[33]:


tested = tested.loc[t,:]


# In[34]:


tested_orfs = tested['orf'].unique()


# In[35]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[36]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[37]:


data = original_data.copy()


# In[38]:


dataset_ids = [16459, 16460]
datasets = datasets.reindex(index=dataset_ids)


# In[39]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[40]:


data.head()


# In[41]:


print((data<0).sum(axis=0))


# ## Subset to the genes currently in SGD

# In[42]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[43]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[44]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[45]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[46]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[47]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[48]:


from IO.save_data_to_db3 import *


# In[49]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




