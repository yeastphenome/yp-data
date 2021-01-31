#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 14718172
paper_name = 'lum_shoemaker_2004' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[25]:


original_data = pd.read_excel('raw_data/mmc2.xlsx', sheet_name='P values')


# In[26]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[27]:


original_data.head()


# In[28]:


original_data['orf'] = original_data['Systematic Name'].astype(str)


# In[29]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[30]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[31]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[32]:


original_data.set_index('orf', inplace=True)


# In[33]:


original_data.drop(columns=['Gene Symbol','Systematic Name'], inplace=True)


# In[34]:


original_data = original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[35]:


original_data = original_data.groupby(original_data.index).mean()


# In[36]:


original_data.shape


# In[50]:


# The entry 'NaN' indicates that no significant growth defect was observed or that data is not available for a particular strain in a particular condition.
# Given the number of NaNs, it is likely that most of them are "no significant growth defect observed", so decided to switch them to 0.
original_data[original_data.isnull()] = 0


# # Load dataset_ids

# In[37]:


dt = pd.read_csv('extras/dataset_ids.txt', sep='\t', header=None)


# In[38]:


dt.head()


# In[39]:


dt.set_index(1, inplace=True)


# In[40]:


dt = dt.reindex(index=original_data.columns.values)


# In[41]:


dt.head()


# In[42]:


dataset_ids = dt[0].values


# # Prepare the final dataset

# In[51]:


data = original_data.copy()


# In[52]:


datasets = datasets.reindex(index=dataset_ids)


# In[53]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[54]:


data.head()


# ## Subset to the genes currently in SGD

# In[55]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[56]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[57]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[58]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[59]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[60]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[61]:


from IO.save_data_to_db3 import *


# In[62]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




