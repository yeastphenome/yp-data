#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 19220866
paper_name = 'mira_sa_correia_2009' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


files = ['hits_genenames.txt','hits_genenames_moderate.txt']


# In[13]:


original_data_list = []
for f in files:
    original_data = pd.read_csv('raw_data/' + f, header=None, sep='\t')
    print('Original data dimensions: %d x %d' % (original_data.shape))
#     print(original_data.head())
    original_data['gene'] = original_data[0].astype(str)
    original_data['gene'] = clean_genename(original_data['gene'])
    original_data['orf'] = translate_sc(original_data['gene'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data = original_data.loc[t,:]
    original_data.set_index('orf', inplace=True)
    original_data['data'] = -1
    original_data = original_data[['data']].copy()
    original_data_list.append(original_data)


# In[14]:


original_data = pd.concat(original_data_list, axis=1)


# In[16]:


original_data.columns=['all','moderate']
original_data[original_data.isnull()] = 0
original_data['data'] = original_data['all'] - original_data['moderate'] - 1


# In[18]:


original_data.index.name='orf'
original_data = original_data[['data']]


# In[19]:


original_data = original_data.groupby(original_data.index).mean()


# In[20]:


original_data.shape


# # Load & process tested strains

# In[21]:


tested = pd.read_excel('raw_data/List of strains tested.xlsx', sheet_name='Tabelle2')


# In[22]:


tested.head()


# In[23]:


tested['orf'] = tested['ORF'].astype(str)


# In[24]:


tested['orf'] = clean_orf(tested['orf'])


# In[25]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[26]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[27]:


tested_orfs = tested['orf'].unique()


# In[28]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[29]:


tested_orfs = list(tested_orfs) + missing


# In[30]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[31]:


data = original_data.copy()


# In[32]:


dataset_ids = [157]
datasets = datasets.reindex(index=dataset_ids)


# In[33]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[34]:


data.head()


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

data.head()


# # Normalize

# In[37]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[38]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[39]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[40]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[41]:


from IO.save_data_to_db3 import *


# In[42]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




