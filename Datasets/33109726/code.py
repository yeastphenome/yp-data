#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[4]:


paper_pmid = 33109726
paper_name = 'ayers_gallagher_2020' 


# In[5]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[6]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[7]:


original_data = pd.read_excel('raw_data/MCHMs KO full list.xlsx', sheet_name='ScreenKOs_withNames_v3')


# In[8]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[9]:


original_data['orfs'] = original_data['Systematic.Name'].astype(str)


# In[10]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[11]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[13]:


original_data = original_data.loc[t,:]


# In[14]:


for c in ['Rep.1','Rep.2','Rep.3','Definite']:
    original_data.loc[original_data.loc[:,c]=='X',c] = -1
    original_data.loc[original_data.loc[:,c].isnull(),c] = 0
    original_data.loc[:,c] = original_data.loc[:,c].astype(int)


# In[15]:


original_data['data'] = original_data[['Rep.1','Rep.2','Rep.3']].sum(axis=1)


# In[16]:


original_data.set_index('orfs', inplace=True)
original_data.index.name='orf'


# In[17]:


original_data.head()


# In[18]:


original_data = original_data.groupby(original_data.index).mean()


# In[19]:


original_data.shape


# In[20]:


original_data.head()


# # Load & process tested strains

# In[21]:


tested = pd.read_excel('raw_data/MATalpha yeast knockout collection.xls', sheet_name='list')


# In[22]:


tested['orfs'] = tested['ORF'].astype(str)


# In[23]:


tested['orfs'] = clean_orf(tested['orfs'])


# In[24]:


tested['orfs'] = translate_sc(tested['orfs'], to='orf')


# In[25]:


# Make sure everything translated ok
t = looks_like_orf(tested['orfs'])
print(tested.loc[~t,])


# In[26]:


tested = tested.loc[t,]


# In[27]:


tested_orfs = np.unique(tested['orfs'].values)


# In[28]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[29]:


tested_orfs = np.append(tested_orfs, 'YPL183W-C')


# In[30]:


# tested_orfs


# In[31]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[32]:


data = original_data[['data']].copy()


# In[33]:


dataset_ids = [16680]
datasets = datasets.reindex(index=dataset_ids)


# In[34]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[35]:


data.head()


# In[36]:


data.shape


# ## Subset to the genes currently in SGD

# In[37]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[38]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[39]:


data.head()


# # Normalize

# In[40]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[41]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[42]:


data_norm[data.isnull()] = np.nan


# In[43]:


data_all = data.join(data_norm)


# In[44]:


data_all.head()


# # Print out

# In[45]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[47]:


from IO.save_data_to_db3 import *


# In[48]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




