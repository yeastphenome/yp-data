#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 26994103
paper_name = 'luo_jiang_2016' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data1 = pd.read_excel('raw_data/Supplementary Table 2.xlsx', sheet_name='Sheet1', skiprows=2, names=['orf','gene','c1','c2'])
original_data2 = pd.read_excel('raw_data/Supplementary Table 3.xlsx', sheet_name='Sheet1', skiprows=2, names=['orf','gene','c1'])


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[7]:


original_data1['orf'] = original_data1['orf'].astype(str)
original_data2['orf'] = original_data2['orf'].astype(str)


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
t1 = looks_like_orf(original_data1['orf'])
t2 = looks_like_orf(original_data2['orf'])


# In[11]:


print(original_data1.loc[~t1,])


# In[12]:


print(original_data2.loc[~t2,])


# In[13]:


original_data1.set_index('orf', inplace=True)
original_data2.set_index('orf', inplace=True)


# In[14]:


original_data1['data'] = -1
original_data2['data'] = 1


# In[15]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[16]:


original_data = original_data[['data_1','data_2']].copy()


# In[17]:


original_data[original_data.isnull()] = 0


# In[18]:


original_data = original_data.groupby(original_data.index).mean()


# In[29]:


original_data.head()


# # Load the tested strains

# In[25]:


tested = pd.read_excel('raw_data/DELETION LIBRARY.xlsx', sheet_name='DELETION LIBRARY', skiprows=1)
tested['orf'] = tested['ORF name'].astype(str)
tested['orf'] = clean_orf(tested['orf'])
tested.loc[tested['orf'] == 'YELOO1C','orf'] = 'YEL001C'
tested['orf'] = translate_sc(tested['orf'], to='orf')

# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])
tested = tested.loc[t,:]


# In[26]:


tested_orfs = np.unique(tested['orf'].values)


# In[27]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[28]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[30]:


data = original_data.copy()


# In[31]:


dataset_ids = [16448,16449]
datasets = datasets.reindex(index=dataset_ids)


# In[32]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[33]:


data.head()


# In[34]:


data.sum(axis=0)


# ## Subset to the genes currently in SGD

# In[35]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[36]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[37]:


data.head()


# # Normalize

# In[38]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[39]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[40]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[41]:


data_all.head()


# # Print out

# In[42]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[43]:


from IO.save_data_to_db3 import *


# In[44]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




