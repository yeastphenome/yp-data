#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 21978764
paper_name = 'svensson_samson_2011' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[28]:


original_data = pd.read_excel('raw_data/1752-0509-5-157-s1.xlsx', sheet_name='2. Gi50 and R2 all strains')


# In[29]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[30]:


original_data.head()


# In[31]:


# Remove the DAMP strains. Personal communication from Peter Svensson: deletions are on plates 1-57, DAMPs are on plates 301-311


# In[32]:


import re


# In[33]:


original_data['plate'] = original_data['position'].apply(lambda x: int(re.findall(r'\d+', x)[0]))


# In[35]:


original_data = original_data.loc[original_data['plate'] <= 57,:]


# In[36]:


original_data['orf'] = original_data['orf'].astype(str)


# In[37]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[38]:


typo_fixes = {'YAR002AW':'YAR002W','YOLO57W':'YOL057W','YKLO72W':'YKL072W',
              'YJL206-A':'YJL206C','YLR287-A':'YLR287C-A','YFL033AC':'YFL033C','YOLO62C':'YOL062C'}
original_data['orf'] = original_data['orf'].apply(lambda x: typo_fixes[x] if x in typo_fixes.keys() else x)


# In[39]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[40]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[41]:


original_data = original_data.loc[t,:]


# In[42]:


original_data['data'] = original_data['GI50']


# In[43]:


original_data.set_index('orf', inplace=True)


# In[44]:


original_data = original_data[['data']].copy()


# In[45]:


original_data = original_data.groupby(original_data.index).mean()


# In[46]:


original_data.shape


# # Prepare the final dataset

# In[47]:


data = original_data.copy()


# In[48]:


dataset_ids = [28]
datasets = datasets.reindex(index=dataset_ids)


# In[49]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[50]:


data.head()


# ## Subset to the genes currently in SGD

# In[51]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[52]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[53]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[54]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[55]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[56]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[57]:


from IO.save_data_to_db3 import *


# In[58]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




