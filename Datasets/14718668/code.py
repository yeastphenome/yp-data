#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 14718668
paper_name = 'giaever_davis_2004' 


# In[103]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[104]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_csv('raw_data/analyzed_data_hom.txt', sep='\t')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data['orf'] = original_data['orf'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[10]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[12]:


original_data.set_index('orf', inplace=True)


# In[14]:


original_data.drop(columns=['gene'], inplace=True)


# In[16]:


original_data = original_data.apply(pd.to_numeric, errors='coerce')


# In[17]:


original_data = original_data.groupby(original_data.index).mean()


# In[18]:


original_data.shape


# # Clean-up

# In[21]:


eliminated = pd.read_csv('raw_data/background_tags_base_hom.txt', sep='\t', header=None)


# In[22]:


eliminated.head()


# In[23]:


eliminated['orf'] = eliminated[0].astype(str)


# In[25]:


eliminated['orf'] = clean_orf(eliminated['orf'])


# In[26]:


eliminated['orf'] = translate_sc(eliminated['orf'], to='orf')


# In[27]:


# Make sure everything translated ok
t = looks_like_orf(eliminated['orf'])
print(eliminated.loc[~t,])


# In[28]:


eliminated_orfs = eliminated['orf'].unique()


# In[29]:


original_data.drop(index=eliminated_orfs, inplace=True)


# In[32]:


# Remove all ORFs that only have 0s because >90% of them are essential genes

n = original_data.sum(axis=1)


# In[35]:


original_data = original_data.loc[n > 0,:]


# In[36]:


original_data.shape


# In[37]:


original_data.columns = [5316, 5316, 5316, 5317, 5317, 5318, 5318, 5313, 5313]


# In[38]:


original_data = original_data.T
original_data = original_data.groupby(original_data.index).mean()
original_data = original_data.T


# # Load HET data

# In[51]:


original_data2 = pd.read_csv('raw_data/analyzed_data_het.txt', sep='\t')


# In[52]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[53]:


original_data2.head()


# In[54]:


original_data2['orf'] = original_data2['orf'].astype(str)


# In[56]:


# Eliminate all white spaces & capitalize
original_data2['orf'] = clean_orf(original_data2['orf'])


# In[57]:


# Translate to ORFs 
original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')


# In[58]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[59]:


original_data2.set_index('orf', inplace=True)


# In[60]:


original_data2.drop(columns=['gene'], inplace=True)


# In[61]:


original_data2 = original_data2.apply(pd.to_numeric, errors='coerce')


# In[62]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[63]:


original_data2.shape


# # Clean-up

# In[64]:


eliminated = pd.read_csv('raw_data/background_tags_base.txt', sep='\t', header=None)


# In[65]:


eliminated.head()


# In[66]:


eliminated['orf'] = eliminated[0].astype(str)


# In[67]:


eliminated['orf'] = clean_orf(eliminated['orf'])


# In[68]:


eliminated['orf'] = translate_sc(eliminated['orf'], to='orf')


# In[69]:


# Make sure everything translated ok
t = looks_like_orf(eliminated['orf'])
print(eliminated.loc[~t,])


# In[70]:


eliminated_orfs = eliminated['orf'].unique()


# In[71]:


original_data2.drop(index=eliminated_orfs, inplace=True)


# In[72]:


het_dataset_ids1 = [5312, 5312, 5312, 5312, 5312, 5314, 5314, 5310, 5315, 5315, 5311, 5337, 5337, 5337, 5337, 5323, 5323, 5322, 5322, 5322]
het_dataset_ids2 = [5321, 5321, 5321, 5319, 5319, 5320, 5320, 5338, 5339, 5339, 5339, 5339]
het_dataset_ids3 = [5334, 5335, 5335, 5336, 5336, 5333, 5333, 5330, 5331, 5331, 5332, 5332, 5325, 5325, 5325, 5324]
het_dataset_ids4 = [484, 484, 484, 5308, 5308, 5308, 5308, 5308, 5308, 5308, 5308, 5308, 5309]
het_dataset_ids5 = [5326, 5326, 5327, 5327, 5327, 5327, 5328, 5328, 5329, 5329]

het_dataset_ids = het_dataset_ids1 + het_dataset_ids2 + het_dataset_ids3 + het_dataset_ids4 + het_dataset_ids5


# In[73]:


original_data2.columns = het_dataset_ids


# In[74]:


original_data2 = original_data2.T
original_data2 = original_data2.groupby(original_data2.index).mean()
original_data2 = original_data2.T


# # Merge

# In[76]:


original_data = original_data.join(original_data2, how='outer')


# In[90]:


# Taking the opposite to make sure that low number indicate low growth and viceversa
original_data = -original_data


# # Prepare the final dataset

# In[105]:


data = original_data.copy()


# In[106]:


dataset_ids = original_data.columns.values
datasets = datasets.reindex(index=dataset_ids)


# In[107]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[108]:


data.head()


# ## Subset to the genes currently in SGD

# In[109]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[110]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[111]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[112]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[113]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[114]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[115]:


from IO.save_data_to_db3 import *


# In[116]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




