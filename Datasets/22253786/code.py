#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 22253786
paper_name = 'blackman_nislow_2012' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/pone.0029798.s003.xlsx', sheet_name = 'elescomol')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data['orf'] = original_data['strain'].apply(lambda x: x.split(':')[0])


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


# In[12]:


original_data['data'] = original_data['log2']


# In[13]:


original_data.set_index('orf', inplace=True)


# In[14]:


original_data = original_data[['data']].copy()


# In[15]:


original_data = original_data.groupby(original_data.index).mean()


# In[16]:


original_data.shape


# # Separate HOM and HET using today's HOM collection from Open Biosystems

# In[19]:


hom = pd.read_excel('extras/Homozygous_diploid_obs_v7.0.xlsx', sheet_name='DATA')


# In[21]:


hom.head()


# In[23]:


hom['orf'] = hom['ORF'].astype(str)


# In[24]:


hom['orf'] = clean_orf(hom['orf'])


# In[25]:


hom['orf'] = translate_sc(hom['orf'], to='orf')


# In[26]:


# Make sure everything translated ok
t = looks_like_orf(hom['orf'])
print(hom.loc[~t,])


# In[27]:


hom_orfs = hom['orf'].unique()


# In[29]:


hom_orfs_in_data = [orf for orf in hom_orfs if orf in original_data.index.values]


# In[31]:


original_data1 = original_data.loc[hom_orfs_in_data,:].copy()
original_data2 = original_data.drop(index=hom_orfs_in_data, inplace=False)


# In[32]:


original_data1.shape


# In[33]:


original_data2.shape


# In[34]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# # Prepare the final dataset

# In[35]:


data = original_data.copy()


# In[36]:


dataset_ids = [86,87]
datasets = datasets.reindex(index=dataset_ids)


# In[37]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[38]:


data.head()


# ## Subset to the genes currently in SGD

# In[39]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[40]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[41]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[42]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[43]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[44]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[45]:


from IO.save_data_to_db3 import *


# In[46]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




