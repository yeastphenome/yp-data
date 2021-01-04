#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 24286259
paper_name = 'sousa_sousa_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[25]:


original_data = pd.read_excel('raw_data/12864_2013_5541_MOESM3_ESM.xlsx', sheet_name='Folha2', skiprows=6)


# In[26]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[27]:


original_data.head(n=10)


# In[28]:


original_data['orf'] = original_data['ORF name'].astype(str)


# In[29]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[39]:


typo_fixes = {'YDL230': 'YDL230W','YDL178': 'YDL178W','YLR124W-': 'YLR124W','YAL016CB': 'YAL016C-B'}
original_data['orf'] = original_data['orf'].apply(lambda x: typo_fixes[x] if x in typo_fixes.keys() else x)


# In[40]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[41]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[43]:


original_data = original_data.loc[t,:]


# In[45]:


data_switch = {'+': 1, '-': -1}
original_data['data'] = original_data['Screening result'].apply(lambda x: data_switch[x])


# In[46]:


original_data.set_index('orf', inplace=True)


# In[47]:


original_data = original_data[['data']].copy()


# In[48]:


original_data = original_data.groupby(original_data.index).mean()


# In[49]:


original_data.shape


# # Prepare the final dataset

# In[50]:


data = original_data.copy()


# In[51]:


dataset_ids = [11812]
datasets = datasets.reindex(index=dataset_ids)


# In[52]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[53]:


data.head()


# ## Subset to the genes currently in SGD

# In[54]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[55]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[56]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[57]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[58]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[59]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[60]:


from IO.save_data_to_db3 import *


# In[61]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




