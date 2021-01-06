#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 23208710
paper_name = 'lis_bobek_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/AAC.01439-12_zac999101582so1.xlsx', sheet_name='all data')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data['orf'] = original_data['ORF'].astype(str)


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


original_data.set_index('orf', inplace=True)


# In[13]:


original_data = -original_data[['MUC7','HSN','KR20','hLF']].copy()


# In[14]:


original_data = original_data.groupby(original_data.index).mean()


# In[15]:


original_data.shape


# In[16]:


original_data.head()


# ### Split essentials and non-essentials

# In[17]:


essential_orfs = original_data.index.values[is_essential(original_data.index.values)]
nonessential_orfs = original_data.index.values[~is_essential(original_data.index.values)]


# In[18]:


original_data1 = original_data.loc[nonessential_orfs,:].copy()
original_data2 = original_data.loc[essential_orfs,:].copy()


# In[19]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[20]:


original_data.head()


# # Prepare the final dataset

# In[21]:


data = original_data.copy()


# In[22]:


dataset_ids = [516, 4925, 4926, 4927, 11870, 11871, 11872, 11873]
datasets = datasets.reindex(index=dataset_ids)


# In[23]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[24]:


data.head()


# ## Subset to the genes currently in SGD

# In[25]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[26]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[27]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[28]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[29]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[30]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[31]:


from IO.save_data_to_db3 import *


# In[32]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




