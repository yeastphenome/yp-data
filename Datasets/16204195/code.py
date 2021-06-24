#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 16204195
paper_name = 'kuepfer_blank_2005' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/15_Juni_APPENDIX.xls', sheet_name='Appendix 1 Mutant library')


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


original_data = original_data.loc[t,:]


# In[13]:


original_data.set_index('orf', inplace=True)


# In[14]:


data_cols = ['glucose(aerobe)', 'glucose (anaerobe)', 'galactose','glycerol','ethanol']
original_data = original_data[data_cols]


# In[15]:


original_data[original_data==0] = -1
original_data[original_data.isnull()] = 0


# In[22]:


for c in original_data.columns:
    original_data[c] = pd.to_numeric(original_data[c], errors='coerce')


# In[23]:


original_data = original_data.groupby(original_data.index).mean()


# In[24]:


original_data.shape


# In[27]:


original_data.sum(axis=0)


# # Prepare the final dataset

# In[28]:


data = original_data.copy()


# In[29]:


dataset_ids = [21954, 21958, 21957, 21955, 21956]
datasets = datasets.reindex(index=dataset_ids)


# In[30]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[31]:


data.head()


# ## Subset to the genes currently in SGD

# In[32]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[33]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[34]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[35]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[36]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[37]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[38]:


from IO.save_data_to_db3 import *


# In[39]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




