#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 27909446
paper_name = 'north_vulpe_2016' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[24]:


original_data = pd.read_excel('raw_data/Table1.xlsx', sheet_name='Strains with altered growth in ', skiprows=10)


# In[25]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[26]:


original_data.head()


# In[27]:


original_data['orf'] = original_data['ORF'].astype(str)


# In[28]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[29]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[30]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[31]:


original_data.set_index('orf', inplace=True)


# In[32]:


original_data = original_data[['150μM','300μM','600μM','150μM.1','300μM.1','600μM.1']].copy()


# In[33]:


original_data = original_data.groupby(original_data.index).mean()


# In[34]:


original_data.shape


# # Prepare the final dataset

# In[35]:


data = original_data.copy()


# In[36]:


dataset_ids = [5377, 5376, 5375, 5380, 5379, 5378]
datasets = datasets.reindex(index=dataset_ids)


# In[37]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[39]:


data[data.isnull()] = 0


# In[40]:


data.head()


# ## Subset to the genes currently in SGD

# In[41]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[42]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[43]:


data.head()


# # Normalize

# In[44]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[45]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[46]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[47]:


data_all.head()


# # Print out

# In[48]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[49]:


from IO.save_data_to_db3 import *


# In[50]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




