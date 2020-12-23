#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 31854076
paper_name = 'novarina_chang_2019' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[10]:


original_data = pd.read_excel('raw_data/acel13084-sup-0001-files1.xlsx', sheet_nam='Filtered', header=8, nrows=4898-9)


# In[11]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[12]:


original_data['orfs'] = original_data['Systematic_Name'].astype(str)


# In[13]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[14]:


original_data.loc[original_data['orfs']=='YLR287-A','orfs'] = 'YLR287C-A'


# In[15]:


original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[16]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[17]:


original_data = original_data.loc[t,:]


# In[18]:


print(original_data.shape)


# In[19]:


original_data.set_index('orfs', inplace=True)
original_data.index.name='orf'


# In[20]:


# Get the relevant columns
original_data = original_data.loc[:, 'percentage_escapers'].to_frame()


# In[21]:


print(original_data.shape)


# In[22]:


original_data = original_data.groupby(original_data.index).mean()


# In[24]:


print('Final data dimensions: %d x %d' % (original_data.shape))


# # Prepare the final dataset

# In[26]:


data = original_data[['percentage_escapers']].copy()


# In[27]:


dataset_ids = [16401]
datasets = datasets.reindex(index=dataset_ids)


# In[28]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[29]:


data.head()


# ## Subset to the genes currently in SGD

# In[30]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[31]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[32]:


data.head()


# # Normalize

# In[33]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[34]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[35]:


data_norm[data.isnull()] = np.nan


# In[36]:


data_all = data.join(data_norm)


# In[37]:


data_all.head()


# # Print out

# In[38]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[39]:


from IO.save_data_to_db3 import *


# In[40]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




