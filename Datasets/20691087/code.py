#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 20691087
paper_name = 'alamgir_golshani_2010' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[26]:


original_data = pd.read_excel('raw_data/1472-6769-10-6-s1.xlsx', sheet_name='Raw genome-wide data', skiprows=1)


# In[27]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[28]:


original_data.head()


# In[29]:


original_data['orf'] = original_data['Systematic Name'].astype(str)


# In[30]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[31]:


original_data.loc[original_data['orf']=='YPL072WA','orf'] = 'YPL072W-A'


# In[32]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[33]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[34]:


original_data = original_data.loc[t,:]


# In[35]:


original_data.set_index('orf', inplace=True)


# In[36]:


original_data = original_data.iloc[:,2:14]


# In[37]:


original_data = -original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[38]:


original_data = original_data.groupby(original_data.index).mean()


# In[39]:


original_data.shape


# In[43]:


original_data.columns = [9, 9, 9, 10, 10, 10, 7, 7, 7, 8, 8, 8]


# In[44]:


original_data = original_data.T
original_data = original_data.groupby(original_data.index).mean()
original_data = original_data.T


# In[45]:


original_data.head()


# # Prepare the final dataset

# In[46]:


data = original_data.copy()


# In[47]:


dataset_ids = [7,8,9,10]
datasets = datasets.reindex(index=dataset_ids)


# In[48]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[49]:


data.head()


# ## Subset to the genes currently in SGD

# In[50]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[51]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[52]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[53]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[54]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[55]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[56]:


from IO.save_data_to_db3 import *


# In[57]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




