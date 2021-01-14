#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 20202201
paper_name = 'batova_schuller_2010' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[53]:


original_data = pd.read_excel('raw_data/1471-2164-11-153-S1.xlsx', sheet_name='Sheet1', header=None)


# In[54]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[55]:


original_data.head()


# In[57]:


t = pd.DataFrame(data = {'orf': original_data.iloc[0::3,0].values, 
                                     'data1': original_data.iloc[1::3,0].values, 
                                     'data2': original_data.iloc[2::3,0].values})


# In[58]:


original_data = t.copy()


# In[59]:


original_data['orf'] = original_data['orf'].astype(str)


# In[60]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[61]:


typo_fixes = {'YJL206CA':'YJL206C-A','YFR031CA':'YFR031C-A','YER087CA':'YER087C-A',
              'YMR031WA':'YMR031W-A','YML081CA':'YML081C-A','YER014CA':'YER014C-A'}


# In[62]:


original_data['orf'] = original_data['orf'].apply(lambda x: typo_fixes[x] if x in typo_fixes.keys() else x)


# In[63]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[64]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[65]:


data_switch = {'+': 0, 'sl': -1, '‐': -2}
for c in ['data1','data2']:
    original_data[c] = original_data[c].apply(lambda x: data_switch[x])


# In[66]:


original_data.head()


# In[67]:


original_data.set_index('orf', inplace=True)


# In[68]:


original_data = original_data.groupby(original_data.index).mean()


# In[69]:


original_data.shape


# # Prepare the final dataset

# In[70]:


data = original_data.copy()


# In[71]:


dataset_ids = [153,437]
datasets = datasets.reindex(index=dataset_ids)


# In[72]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[73]:


data.head()


# ## Subset to the genes currently in SGD

# In[74]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[75]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[79]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[80]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[81]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[82]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[83]:


from IO.save_data_to_db3 import *


# In[84]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




