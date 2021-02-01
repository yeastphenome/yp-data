#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 14594803
paper_name = 'davis_kaplan_kaplan_2004' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/Copy of Kaplan Experimental diploid homozygous deletion srceen Release 1 and 2 sort sheet(3).xlsx', sheet_name='Exper. Results 1 & 2 sort sheet')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data['orf'] = original_data[' ORF'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[10]:


typo_fixes = pd.read_csv('extras/typo_fixes.txt', sep='\t', header=None)
typo_fixes_dict = {row[0]: row[1] for index, row in typo_fixes.iterrows() }
original_data['orf'] = original_data['orf'].apply(lambda x: typo_fixes_dict[x] if x in typo_fixes_dict.keys() else x)


# In[11]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[13]:


original_data = original_data.loc[t,:]


# In[14]:


original_data.set_index('orf', inplace=True)


# In[15]:


data_cols = [c for c in original_data.columns.values if c.startswith('3')]


# In[16]:


original_data = original_data[data_cols].copy()


# In[17]:


original_data = original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[18]:


original_data = original_data.groupby(original_data.index).mean()


# In[19]:


original_data.shape


# In[20]:


cols = [c.strip('.1').replace('°','C').replace('µ','u').replace(' ','') for c in original_data.columns.values]


# # Load dataset ids

# In[21]:


dt = pd.read_csv('extras/dataset_ids.txt', sep='\t', header=None)


# In[22]:


dt.head()


# In[23]:


dt.shape


# In[24]:


dt[0] = dt[0].apply(lambda x: x.replace('°','C').replace('µ','u').replace(' ',''))


# In[25]:


dt.set_index(0, inplace=True)


# In[26]:


dt.head()


# In[27]:


dataset_ids = [dt.loc[x,1] for x in cols]


# In[28]:


len(dataset_ids)


# In[29]:


original_data.columns = dataset_ids


# In[30]:


original_data = original_data.T
original_data = original_data.groupby(original_data.index).mean()
original_data = original_data.T


# In[31]:


## Exclude the datasets with < 1000 tested genes
num_tested = np.sum(original_data.notnull(), axis=0)


# In[32]:


original_data = original_data.loc[:, num_tested >= 1000]


# In[33]:


original_data.shape


# # Prepare the final dataset

# In[34]:


data = original_data.copy()


# In[35]:


dataset_ids = original_data.columns.values
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


# In[39]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

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
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[43]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[44]:


from IO.save_data_to_db3 import *


# In[45]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




