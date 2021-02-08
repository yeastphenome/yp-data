#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[34]:


paper_pmid = 12543677
paper_name = 'blackburn_avery_2003' 


# In[35]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[36]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[37]:


original_data = pd.read_excel('raw_data/blackburn_avery_2003_data.xlsx', sheet_name='data.txt')


# In[38]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[39]:


original_data.head()


# In[40]:


original_data['orf'] = original_data['Unnamed: 0'].astype(str)


# In[41]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[42]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[43]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[44]:


original_data.set_index('orf', inplace=True)


# In[45]:


# Removing 4 datasets because they don't didn't produce any hits
original_data.drop(columns=['Unnamed: 0','amoxicillin','penicillin G','rifampin','vancomycin'], inplace=True)


# In[46]:


original_data = original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[47]:


# MIC=Inf set to 600 (+88 relative to the maximum MIC detected)
vals = original_data.values
vals[np.isinf(vals)] = 600
vals = vals - 600 # Transform so that lack of sensivity (=wt) = 0
original_data = pd.DataFrame(index=original_data.index, columns=original_data.columns, data=vals)


# In[48]:


original_data = original_data.groupby(original_data.index).mean()


# In[49]:


original_data.shape


# In[50]:


np.abs(original_data).sum(axis=0)


# # Prepare the final dataset

# In[51]:


data = original_data.copy()


# In[52]:


dataset_ids = [391, 393, 395]
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


data_norm = normalize_phenotypic_scores(data, has_tested=False)


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




