#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 12477387
paper_name = 'zhang_schneider_2002' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[8]:


original_data = pd.read_excel('raw_data/Zhang et al supplemental data.xlsx', sheet_name='Size Data', header=None)


# In[9]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[10]:


original_data.head()


# In[11]:


original_data['orf'] = original_data[0].astype(str)


# In[12]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[15]:


typo_fixes = {'TAL004W':'YAL004W','YELOO1C':'YEL001C','KL187C':'YKL187C'}
original_data['orf'] = original_data['orf'].apply(lambda x: typo_fixes[x] if x in typo_fixes.keys() else x)


# In[16]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[17]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[18]:


original_data = original_data.loc[t,:]


# In[19]:


original_data['data'] = pd.to_numeric(original_data[4], errors='coerce')


# In[20]:


original_data.set_index('orf', inplace=True)


# In[21]:


original_data = original_data[['data']].copy()


# In[22]:


original_data = original_data.groupby(original_data.index).mean()


# In[23]:


original_data.shape


# In[24]:


# To separate HOM from HET data, use genes on today's HOM collection from Open Biosystems
hom = pd.read_excel('extras/Homozygous_diploid_obs_v7.0.xlsx', sheet_name='DATA')
hom.head()


# In[25]:


hom['orf'] = hom['ORF'].astype(str)


# In[26]:


hom['orf'] = clean_orf(hom['orf'])


# In[27]:


hom['orf'] = translate_sc(hom['orf'], to='orf')


# In[28]:


t = looks_like_orf(hom['orf'])
print(hom.loc[~t,])


# In[29]:


hom_orfs = hom['orf'].unique()


# In[30]:


orfs_in_hom = [orf for orf in original_data.index.values if orf in hom_orfs]
orfs_in_het = [orf for orf in original_data.index.values if orf not in hom_orfs]


# In[31]:


original_data1 = original_data.loc[orfs_in_hom,:].copy()
original_data2 = original_data.loc[orfs_in_het,:].copy()


# In[32]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[33]:


original_data.head()


# # Prepare the final dataset

# In[34]:


data = original_data.copy()


# In[35]:


dataset_ids = [477, 5384]
datasets = datasets.reindex(index=dataset_ids)


# In[36]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[37]:


data.head()


# ## Subset to the genes currently in SGD

# In[38]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[41]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[42]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[43]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[44]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[45]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[46]:


from IO.save_data_to_db3 import *


# In[47]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




