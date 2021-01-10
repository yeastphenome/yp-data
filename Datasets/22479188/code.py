#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 22479188
paper_name = 'albulescu_pleiss_2012' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[29]:


original_data = pd.read_excel('raw_data/hits.xlsx', sheet_name='Sheet1')


# In[30]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[31]:


original_data.head()


# In[32]:


original_data['orf'] = original_data['NAME'].astype(str)


# In[33]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[34]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[35]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[36]:


original_data.set_index('orf', inplace=True)


# In[37]:


original_data = original_data[['U3','Tef5','Rpl31b','Tub3','Ubc13']].copy()


# In[38]:


original_data = original_data.groupby(original_data.index).mean()


# In[39]:


original_data.shape


# In[41]:


original_data[original_data.isnull()] = 1


# In[42]:


original_data.head()


# # Load & process tested strains

# In[43]:


tested = pd.read_csv('raw_data/no_hits.txt', sep='\t')


# In[44]:


tested.head()


# In[45]:


tested['orf'] = tested['Systematic name'].astype(str)


# In[46]:


tested['orf'] = clean_orf(tested['orf'])


# In[47]:


tested.loc[tested['orf']=='YLR287-A','orf'] = 'YLR287C-A'
tested.loc[tested['orf']=='YCRO54C','orf'] = 'YCR054C'


# In[48]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[49]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[50]:


tested = tested.loc[t,:]


# In[51]:


tested_orfs = tested['orf'].unique()


# In[52]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[53]:


original_data = original_data.reindex(index=tested_orfs, fill_value=1)


# # Prepare the final dataset

# In[54]:


data = original_data.copy()


# In[55]:


dataset_ids = [16311, 16315, 16312, 16314, 16313]
datasets = datasets.reindex(index=dataset_ids)


# In[56]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[57]:


data.head()


# ## Subset to the genes currently in SGD

# In[58]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[59]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[60]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[61]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[62]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[63]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[64]:


from IO.save_data_to_db3 import *


# In[65]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




