#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 32919014
paper_name = 'guan_zhang_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/Table S1 Normolized counts.xlsx', sheet_name='69samples normolized bigger tha', skiprows=1)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data['orf'] = original_data['Unnamed: 0'].astype(str)


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


# In[ ]:





# In[16]:


original_data.drop(columns=['Unnamed: 0'], inplace=True)


# In[12]:


original_data.set_index('orf', inplace=True)


# In[17]:


cols = ['_'.join(c.split("_")[0:-1]) for c in original_data.columns]


# In[23]:


original_data.columns = cols


# In[24]:


original_data = original_data.T


# In[25]:


original_data = original_data.groupby(original_data.index.values).mean().T


# In[26]:


original_data.head()


# In[27]:


original_data.shape


# In[31]:


orfs_essential = original_data.index.values[is_essential(original_data.index.values)]
orfs_nonessential = original_data.index.values[~is_essential(original_data.index.values)]


# In[32]:


original_data1 = original_data.loc[orfs_nonessential,:].copy()
original_data2 = original_data.loc[orfs_essential,:].copy()


# In[39]:


for c in original_data1.columns:
    if not 'DMSO' in c:
        original_data1[c] = original_data1[c] / original_data1['DMSO']


# In[40]:


for c in original_data2.columns:
    if not 'DMSO' in c:
        original_data2[c] = original_data2[c] / original_data2['DMSO']


# In[54]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_noness', rsuffix='_ess')


# In[55]:


original_data.shape


# In[56]:


original_data.head()


# # Prepare the final dataset

# In[57]:


data = original_data.copy()


# In[58]:


original_data.columns


# In[59]:


dataset_ids_noness = [21931, 21886, 21921, 21922, 21885, 21927, 21928, 21883, 21923, 21924, 21882, 21925, 21926, 21884, 21929, 21930]
dataset_ids_ess = [21934, 21937, 21936, 21935, 21940, 21938, 21939, 21943, 21942, 21941, 21946, 21945, 21944, 21949, 21948, 21947]
dataset_ids = dataset_ids_noness + dataset_ids_ess
datasets = datasets.reindex(index=dataset_ids)


# In[60]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[61]:


data.head()


# ## Subset to the genes currently in SGD

# In[62]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[63]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# In[64]:


data.shape


# # Normalize

# In[65]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[66]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[67]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[68]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[69]:


from IO.save_data_to_db3 import *


# In[70]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




