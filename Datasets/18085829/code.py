#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 18085829
paper_name = 'alvaro_rothstein_2007' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[14]:


original_data = pd.read_excel('raw_data/pgen.0030228.st001.xlsx', sheet_name='SIA-Rad52 focus screen', skiprows=7)


# In[15]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[16]:


original_data.head()


# In[17]:


original_data['orf'] = original_data['ORF'].astype(str)


# In[18]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[21]:


original_data.loc[original_data['Alias']=='YLR006C','orf'] = 'YLR006C'


# In[22]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[23]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[24]:


original_data['data'] = pd.to_numeric(original_data['% Cells with foci'], errors='coerce')


# In[25]:


original_data.set_index('orf', inplace=True)


# In[26]:


original_data = original_data[['data']].copy()


# In[27]:


original_data = original_data.groupby(original_data.index).mean()


# In[28]:


original_data.shape


# # Prepare the final dataset

# In[29]:


data = original_data.copy()


# In[30]:


dataset_ids = [591]
datasets = datasets.reindex(index=dataset_ids)


# In[31]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[32]:


data.head()


# ## Subset to the genes currently in SGD

# In[33]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[34]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[35]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[36]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[37]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

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



