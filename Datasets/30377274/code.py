#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 30377274
paper_name = 'perez_samper_verstrepen_2018' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/RawData_BarSeq_Gemma.xlsx', sheet_name='Calculate_extendedGRs', skiprows=23)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[17]:


original_data['orf'] = original_data['NameGen'].astype(str)


# In[18]:


original_data['orf'] = original_data['orf'].apply(lambda x: x.replace('_','-'))


# In[19]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[20]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[21]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[22]:


data_cols = [c for c in original_data.columns if c.startswith('GR')]


# In[23]:


data_cols


# In[24]:


original_data.set_index('orf', inplace=True)


# In[25]:


original_data = original_data[data_cols].copy()


# In[29]:


# Average A and B, Up and Down
cs = ['_'.join(c.split('_')[:3]) for c in data_cols]


# In[31]:


original_data.columns = cs
original_data = original_data.groupby(original_data.columns, axis=1).mean()


# In[32]:


original_data.shape


# In[33]:


original_data = original_data.groupby(original_data.index).mean()


# In[34]:


original_data.shape


# In[35]:


original_data.head()


# In[48]:


# Normalize gradual shift to glucose 5%
original_data['GR_Grad_Gal'] = original_data['GR_Grad_Gal'] / original_data['GR_Stab_Glu']


# # Prepare the final dataset

# In[49]:


data = original_data.copy()


# In[50]:


dataset_ids = [21907, 21961, 21960]
datasets = datasets.reindex(index=dataset_ids)


# In[51]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[52]:


data.head()


# ## Subset to the genes currently in SGD

# In[53]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[54]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[55]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[56]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[57]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[58]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[59]:


from IO.save_data_to_db3 import *


# In[60]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




