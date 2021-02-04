#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 12134146
paper_name = 'steinmetz_davis_2002' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[41]:


files = ['Regression_Tc1_hom.txt','Regression_Tc2_hom.txt','Regression_Tc1_het.txt','Regression_Tc2_het.txt']


# In[42]:


original_data_list = []
for f in files:
    original_data = pd.read_csv('raw_data/' + f, sep='\t')
    print('Original data dimensions: %d x %d' % (original_data.shape))
#     print(original_data1.head())
    original_data['orf'] = original_data['orf'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data.set_index('orf', inplace=True)
    
    data_cols = original_data.columns[6:]
    original_data = original_data[data_cols].apply(pd.to_numeric, axis=1, errors='coerce')
    original_data = original_data.groupby(original_data.index).mean()
    
    suffix = '_het' if f.endswith('het.txt') else '_hom'
    cols = [c+suffix for c in original_data.columns]
    original_data.columns = cols
    
    original_data_list.append(original_data)


# In[43]:


original_data = pd.concat(original_data_list, axis=1)


# In[44]:


original_data.head()


# In[45]:


original_data = original_data.T
original_data = original_data.groupby(original_data.index).mean()
original_data = original_data.T


# In[46]:


original_data.head()


# In[50]:


# Remove the rows that are all NaNs
num_vals = original_data.notnull().sum(axis=1)


# In[48]:


num_vals.unique()


# In[51]:


original_data = original_data.loc[num_vals > 0,:]


# In[52]:


original_data.shape


# In[53]:


# Normalize all datasets to YPD
# SAFE analysis indicates that taking the difference makes more sense than the ratio
cols_hom = [c for c in original_data.columns if c.endswith('_hom')]
original_data[cols_hom] = original_data[cols_hom].subtract(original_data['YPD_hom'], axis=0)


# In[54]:


cols_het = [c for c in original_data.columns if c.endswith('_het')]
original_data[cols_het] = original_data[cols_het].subtract(original_data['YPD_het'], axis=0)


# In[55]:


original_data.drop(columns=['YPD_hom','YPD_het'], inplace=True)


# In[56]:


original_data.shape


# In[59]:


dataset_ids = [4831, 4841, 4832, 4842, 4830, 4840, 4828, 4838, 4833, 4843, 4834, 4844, 4827, 4837, 4829, 4839]


# In[60]:


original_data.columns = dataset_ids


# # Prepare the final dataset

# In[62]:


data = original_data.copy()


# In[63]:


dataset_ids = original_data.columns.values
datasets = datasets.reindex(index=dataset_ids)


# In[64]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[65]:


data.head()


# ## Subset to the genes currently in SGD

# In[66]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[67]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[68]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[69]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[70]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[71]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[72]:


from IO.save_data_to_db3 import *


# In[73]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




