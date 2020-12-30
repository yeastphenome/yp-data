#!/usr/bin/env python
# coding: utf-8

# In[24]:


get_ipython().run_line_magic('run', '../yp_utils.py')
import itertools


# # Initial setup

# In[2]:


paper_pmid = 26324718
paper_name = 'deranieh_greenberg_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[8]:


original_data = pd.read_csv('raw_data/hit_list.txt', sep='\t', header=None)


# In[9]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[22]:


original_data[3] = original_data[2].apply(lambda x: [g.strip().split('/')[0] for g in x.split(',')])


# In[23]:


original_data.head()


# In[26]:


hit_list = list(itertools.chain.from_iterable(original_data[3]))


# In[30]:


original_data2 = pd.DataFrame(data={'genes': hit_list,'data': np.zeros(len(hit_list))-1})


# In[32]:


original_data2['genes'] = original_data2['genes'].astype(str)


# In[33]:


# Eliminate all white spaces & capitalize
original_data2['genes'] = clean_genename(original_data2['genes'])


# In[34]:


# Translate to ORFs 
original_data2['orf'] = translate_sc(original_data2['genes'], to='orf')


# In[35]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[36]:


original_data2 = original_data2.loc[t,:]


# In[37]:


original_data2.set_index('orf', inplace=True)


# In[38]:


original_data2 = original_data2[['data']].copy()


# In[39]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[40]:


original_data2.shape


# # Prepare the final dataset

# In[41]:


data = original_data2.copy()


# In[42]:


dataset_ids = [4955]
datasets = datasets.reindex(index=dataset_ids)


# In[43]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[44]:


data.head()


# ## Subset to the genes currently in SGD

# In[45]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[46]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[47]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[48]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[49]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[50]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[51]:


from IO.save_data_to_db3 import *


# In[52]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




