#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 29329531
paper_name = 'bottoms_piotrowski_2018' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[20]:


original_data = pd.read_excel('raw_data/YPDg+2.3%gval.v.YPD+gval.pw.rc.edger.2013_12_01_2 for spotfire 001.xlsx', 
                              sheet_name='YPDg+2.3%gval.v.YPD+gval.pw.rc.')


# In[21]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[22]:


original_data.head()


# In[23]:


original_data['genes'] = original_data['ORF'].astype(str)


# In[24]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[25]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[26]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[27]:


# Untranslated ORFs are all pseudogenes or blocked ORFs
original_data = original_data.loc[t,:]


# In[28]:


original_data['data'] = original_data['foldChange']


# In[29]:


original_data.set_index('orfs', inplace=True)
original_data.index.name='orf'


# In[31]:


original_data = original_data.groupby(original_data.index).mean()


# In[32]:


original_data.shape


# In[33]:


# This dataset (for some reason) contains 573 essential genes -- to remove


# In[34]:


essential_orfs = original_data.index.values[is_essential(original_data.index.values)]


# In[35]:


original_data.drop(index=essential_orfs, inplace=True)


# In[37]:


original_data.shape


# # Prepare the final dataset

# In[38]:


data = original_data[['data']].copy()


# In[39]:


dataset_ids = [16211]
datasets = datasets.reindex(index=dataset_ids)


# In[40]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[41]:


data.head()


# ## Subset to the genes currently in SGD

# In[42]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[43]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[44]:


data.head()


# # Normalize

# In[45]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[46]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[47]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[48]:


data_all.head()


# # Print out

# In[49]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[50]:


from IO.save_data_to_db3 import *


# In[51]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




