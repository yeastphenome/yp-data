#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 26240381
paper_name = 'tigano_ottonello_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[8]:


original_data = pd.read_excel('raw_data/nar-03484-z-2014-File007.xlsx', sheet_name='Suppl_TabS1', skiprows=2)


# In[9]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[10]:


original_data.head()


# In[11]:


original_data['orf'] = original_data['Systematic name'].astype(str)


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


original_data.set_index('orf', inplace=True)


# In[17]:


data_switch = {'HS': -3, 'MS': -2, 'LS': -1}
data_cols = ['Ethanol, 30°C','Ethanol, 37°C']

for c in data_cols:
    original_data[c] = original_data[c].apply(lambda x: data_switch[x] if x in data_switch.keys() else x)


# In[18]:


original_data = original_data[data_cols].copy()


# In[20]:


original_data = original_data.groupby(original_data.index).mean()


# In[21]:


original_data.shape


# # Load & process tested strains

# In[22]:


tested = pd.read_excel('raw_data/Tigano et al NAR 2015_List of tested yeast mutant strains.xlsx', sheet_name='Foglio1')


# In[23]:


tested.head()


# In[24]:


tested['orf'] = tested['ORF name'].astype(str)
tested['orf'] = clean_orf(tested['orf'])


# In[25]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[26]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[27]:


tested_orfs = np.unique(tested['orf'].values)


# In[28]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[29]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[30]:


data = original_data.copy()


# In[31]:


dataset_ids = [695, 694]
datasets = datasets.reindex(index=dataset_ids)


# In[32]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[33]:


data.head()


# ## Subset to the genes currently in SGD

# In[34]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[35]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[36]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[37]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[38]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[39]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[40]:


from IO.save_data_to_db3 import *


# In[41]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:



