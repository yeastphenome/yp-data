#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 16582425
paper_name = 'hancock_lopes_2006' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_csv('raw_data/hits_genenames.txt', header=None, names=['gene'])


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data['gene'] = original_data['gene'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['gene'] = clean_genename(original_data['gene'])


# In[10]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['gene'], to='orf')


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[12]:


original_data['data'] = 1


# In[13]:


original_data.set_index('orf', inplace=True)


# In[14]:


original_data = original_data[['data']].copy()


# In[15]:


original_data = original_data.groupby(original_data.index).mean()


# In[16]:


original_data.shape


# # Load & process tested strains

# In[19]:


tested = pd.read_excel('raw_data/mat_alpha_061101.xlsx', sheet_name='mat_alpha_061101', skiprows=1)


# In[20]:


tested.head()


# In[21]:


tested['orf'] = tested['ORF name'].astype(str)


# In[22]:


tested['orf'] = clean_orf(tested['orf'])


# In[25]:


tested.loc[tested['orf']=='YYKL138C','orf'] = 'YKL138C'


# In[26]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[27]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[28]:


tested = tested.loc[t,:]


# In[29]:


tested_orfs = tested['orf'].unique()


# In[30]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[31]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[32]:


data = original_data.copy()


# In[33]:


dataset_ids = [178]
datasets = datasets.reindex(index=dataset_ids)


# In[34]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[35]:


data.head()


# ## Subset to the genes currently in SGD

# In[36]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[37]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

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

data_all.head()


# # Print out

# In[41]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[42]:


from IO.save_data_to_db3 import *


# In[43]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:



