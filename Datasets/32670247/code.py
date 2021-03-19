#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[5]:


paper_pmid = 32670247
paper_name = 'johnston_strobel_2020' 


# In[6]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[7]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[8]:


original_data = pd.read_excel('raw_data/Table_S1.xlsx', sheet_name='Sheet1')


# In[9]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[10]:


original_data.head()


# In[11]:


original_data['orf'] = original_data['ORF name'].astype(str)


# In[12]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[13]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[14]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[15]:


original_data = original_data.loc[t,:]


# In[16]:


original_data.set_index('orf', inplace=True)


# In[17]:


data_cols = ['FCCP','DNP','H2SO4','HCl','NaF']


# In[18]:


original_data = original_data[data_cols].copy()


# In[22]:


data_dict = {'X': -2, 'x': -1, np.nan: 0}
for c in data_cols:
    original_data[c] = original_data[c].apply(lambda x: data_dict[x])


# In[23]:


original_data.head()


# In[24]:


original_data = original_data.groupby(original_data.index).mean()


# In[25]:


original_data.shape


# # Load & process tested strains

# In[26]:


tested = pd.read_excel('raw_data/Matav50.xlsx', sheet_name='DATA')


# In[27]:


tested.head()


# In[33]:


tested['orf'] = tested['ORF name'].astype(str)


# In[34]:


tested['orf'] = clean_orf(tested['orf'])


# In[35]:


tested.loc[tested['orf']=='YLR287-A','orf'] = 'YLR287C-A'


# In[36]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[37]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[38]:


tested = tested.loc[t,:]


# In[39]:


tested_orfs = tested['orf'].unique()


# In[40]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[41]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# In[43]:


original_data.head()


# In[44]:


# Split data into 2 screens per treatment: LOAEL (lower dose) and IC25 (higher dose)
original_data[['FCCP2','DNP2','H2SO42','HCl2','NaF2']] = original_data[['FCCP','DNP','H2SO4','HCl','NaF']]


# In[51]:


low_dose_cols = ['FCCP','DNP','H2SO4','HCl','NaF']
high_dose_cols = ['FCCP2','DNP2','H2SO42','HCl2','NaF2']

t1 = original_data[low_dose_cols].mask(original_data[low_dose_cols]<0, other=-1)
t2 = original_data[high_dose_cols].mask(original_data[high_dose_cols]>-2, other=0)
t2 = t2.mask(t2<0, other=-1)


# In[52]:


original_data = t1.join(t2, how='outer')


# # Prepare the final dataset

# In[58]:


data = original_data.copy()


# In[59]:


dataset_ids = [21853, 21855, 21861, 21859, 21857, 21852, 21854, 21860, 21858, 21856]
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


# # Normalize

# In[64]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[65]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[66]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[67]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[68]:


from IO.save_data_to_db3 import *


# In[69]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




