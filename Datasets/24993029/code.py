#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 24993029
paper_name = 'walker_jiranek_2014' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[77]:


original_data = pd.read_excel('raw_data/12864_2013_6243_MOESM1_ESM.xlsx', sheet_name="Add' file 1 BMC Genomics 2013", skiprows=10)


# In[78]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[79]:


original_data.head()


# In[80]:


original_data['orf'] = original_data['Systematic Name'].astype(str)


# In[81]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[82]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[83]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[84]:


original_data['data'] = pd.to_numeric(original_data['Average'], errors='coerce')


# In[85]:


original_data.loc[original_data['orf']=='BY4743','data']


# In[90]:


# Normalize to WT
original_data['data'] = original_data.loc[original_data['orf']=='BY4743','data'].values - original_data['data'].values


# In[92]:


original_data = original_data.loc[t,:]


# In[93]:


original_data.set_index('orf', inplace=True)


# In[94]:


original_data = original_data[['data']].copy()


# In[95]:


original_data = original_data.groupby(original_data.index).mean()


# In[96]:


original_data.shape


# # Prepare the final dataset

# In[97]:


data = original_data.copy()


# In[98]:


dataset_ids = [11811]
datasets = datasets.reindex(index=dataset_ids)


# In[99]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[100]:


data.head()


# ## Subset to the genes currently in SGD

# In[101]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[102]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[103]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[104]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[105]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[106]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[107]:


from IO.save_data_to_db3 import *


# In[108]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




