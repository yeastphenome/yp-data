#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 31029968
paper_name = 'alfatah_arumugam_2019' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[24]:


original_data = pd.read_excel('raw_data/Supplementary table. 1.xlsx', sheet_name='HOP_Monoethylhexylphthalicacid_')


# In[25]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[26]:


original_data['genes'] = original_data['genes'].astype(str)


# In[27]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[28]:


# If possible, fix typos, omissions, etc.
original_data.loc[original_data['genes'].str.contains('2001-10-01'),'genes'] = 'OCT1'
original_data.loc[original_data['genes'].str.contains('YBR160WAS'),'genes'] = 'YBR160W'


# In[29]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[30]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])


# In[31]:


print(original_data.loc[~t,])


# In[32]:


# If the same strain is present more than once, average its values
data = original_data.groupby('orfs')['logFC'].mean().to_frame()


# In[33]:


print('Final data dimensions: %d x %d' % (data.shape))


# In[34]:


data.index.name='orf'


# # Prepare the final dataset

# In[35]:


dataset_ids = [16439]
datasets = datasets.reindex(index=dataset_ids)


# In[36]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[37]:


data.shape


# In[38]:


data.head()


# ## Subset to the genes currently in SGD

# In[39]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[40]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[41]:


data.head()


# # Normalize

# In[42]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[43]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[44]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[45]:


data_all.head()


# # Print out

# In[46]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[47]:


from IO.save_data_to_db3 import *


# In[48]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




