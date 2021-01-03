#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 24926745
paper_name = 'tun_wu_2014' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[13]:


original_data = pd.read_excel('raw_data/c4mt00116h1.xlsx', sheet_name='Sheet1', skiprows=1)


# In[14]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[15]:


original_data.head()


# In[16]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[17]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[18]:


original_data.loc[original_data['ORF'].str.startswith('YOR205CHOMDIP'),'ORF'] = 'YOR205C'


# In[19]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'], to='orf')


# In[20]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[22]:


original_data = original_data[['ORF','Control','Al 1.6 mM','Al  3.2 mM']].copy()


# In[24]:


original_data.set_index('ORF', inplace=True)
original_data.index.name = 'orf'


# In[25]:


original_data['Control'] = pd.to_numeric(data['Control'], errors='coerce')
original_data['Al 1.6 mM'] = pd.to_numeric(data['Al 1.6 mM'], errors='coerce')
original_data['Al  3.2 mM'] = pd.to_numeric(data['Al  3.2 mM'], errors='coerce')


# In[26]:


original_data = original_data.div(original_data.loc['BY4743',:])


# In[27]:


original_data['Al 1.6 mM'] = original_data['Al 1.6 mM'] / original_data['Control']


# In[28]:


original_data['Al  3.2 mM'] = original_data['Al  3.2 mM'] / original_data['Control']


# In[29]:


original_data.drop(index='BY4743', inplace=True)


# In[30]:


original_data.head()


# In[31]:


original_data = original_data.groupby(original_data.index).mean()


# In[32]:


original_data.shape


# # Prepare the final dataset

# In[33]:


data = original_data.copy()


# In[34]:


dataset_ids = [16509,16477,16478]
datasets = datasets.reindex(index=dataset_ids)


# In[35]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[36]:


data.head()


# ## Subset to the genes currently in SGD

# In[37]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[38]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[39]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[40]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[41]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[42]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[43]:


from IO.save_data_to_db3 import *


# In[44]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




