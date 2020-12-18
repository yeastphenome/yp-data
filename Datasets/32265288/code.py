#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')

import regex


# # Initial setup

# In[2]:


paper_pmid = 32265288
paper_name = 'novarina_chang_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data1 = pd.read_excel('raw_data/Supp_Tables_v4.xlsx', sheet_name='Table S2', skiprows=2)
original_data2 = pd.read_excel('raw_data/Supp_Tables_v4.xlsx', sheet_name='Table S4', skiprows=2)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))


# In[7]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# ### Dataset 1

# In[8]:


gene_col1 = 'Positives from patch assay'


# In[9]:


original_data1['genes'] = original_data1[gene_col1].astype(str)


# In[10]:


# Eliminate all white spaces & capitalize
original_data1['genes'] = clean_genename(original_data1['genes'])


# In[11]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['genes'], to='orf')


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[13]:


original_data1 = original_data1.loc[t,:]


# In[14]:


original_data1['data'] = 1


# In[15]:


original_data1.set_index('orfs', inplace=True)


# In[16]:


original_data1 = original_data1.groupby(original_data1.index).mean()


# ### Dataset 2

# In[17]:


orf_col = 'ORF name'


# In[18]:


original_data2['orfs'] = original_data2[orf_col].astype(str)


# In[19]:


# Eliminate all white spaces & capitalize
original_data2['orfs'] = clean_orf(original_data2['orfs'])


# In[20]:


# Translate to ORFs 
original_data2['orfs'] = translate_sc(original_data2['orfs'], to='orf')


# In[21]:


original_data2.loc[original_data2['orfs']=='YLR287-A','orfs'] = 'YLR287C-A'


# In[22]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[23]:


wt_phenotype = 56 # from main text
original_data2['data'] = (original_data2['Percent recombinants'] / wt_phenotype) - 1


# In[24]:


original_data2.set_index('orfs', inplace=True)


# In[25]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# ### Join the 2 datasets

# In[26]:


original_data = original_data1[['data']].join(original_data2[['data']], how='outer', lsuffix='_1', rsuffix='_2')


# In[27]:


original_data.shape


# In[28]:


# Orfs from dataset 1 that were not present in dataset 2
unmatched_orfs = original_data.loc[original_data['data_2'].isnull()].index.values


# In[29]:


# Check for partial matches of hits from dataset1 in the list from dataset2 (possible typos in dataset1)
orfs_ref = pd.Series(original_data2.index.values)

for orf in unmatched_orfs:
    s = '(' + orf + '){e<=1}'
    
    m = orfs_ref.apply(lambda x: len(regex.findall(s, x)))
    nm = m.sum()
    
    print('%s\t%d' % (orf, nm))


# In[30]:


# No obvious typo fixes. Decided to leave dataset 2 values for unmatched orfs at NaN. But use dataset 2 list as the tested universe for dataset 1 (best approximation).
original_data['data_1'].loc[original_data['data_1'].isnull()] = 0


# In[31]:


original_data.notnull().sum(axis=0)


# # Prepare the final dataset

# In[32]:


dataset_ids = [16617, 16618]


# In[33]:


datasets = datasets.reindex(index=dataset_ids)


# In[34]:


data = original_data[['data_1','data_2']].copy()


# In[35]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[36]:


data = data.groupby(data.index).mean()


# In[37]:


# Create row index
data.index.name='orf'


# In[38]:


print('Final data dimensions: %d x %d' % (data.shape))


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


# In[44]:


data_all = data.join(data_norm)


# In[45]:


data_all.head()


# # Print out

# In[46]:


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




