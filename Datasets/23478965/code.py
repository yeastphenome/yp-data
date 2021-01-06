#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 23478965
paper_name = 'richie_hoepfner_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[29]:


files = ['HOP-exp-scores-annotation.txt','HIP-exp-scores-annotation.txt']


# In[30]:


original_data_list = []
for f in files:
    original_data = pd.read_csv('large_files/raw_data/' + f, sep='\t')
    print('Original data dimensions: %d x %d' % (original_data.shape))
    original_data = pd.pivot_table(original_data, index='Systematic Names', columns='Compound_conc', values='Score', fill_value=np.nan)
    print(original_data.shape)
    
    original_data['orf'] = original_data.index.values
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    original_data.loc['YPR099C','orf'] = 'YPR099C'
    
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data = original_data.loc[t,:]
    original_data.set_index('orf', inplace=True)
    original_data = original_data.groupby(original_data.index).mean()
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[31]:


original_data = original_data_list[0].join(original_data_list[1], how='outer', lsuffix='_1', rsuffix='_2')


# In[32]:


original_data.head()


# # Prepare the final dataset

# In[33]:


data = original_data.copy()


# In[34]:


dataset_ids_hop = [16018, 16008, 16009, 16010, 16011, 16012, 16013, 16014, 16015, 16016, 16017]
dataset_ids_hip = [16029, 16019, 16020, 16021, 16022, 16023, 16024, 16025, 16026, 16027, 16028]
dataset_ids = dataset_ids_hop + dataset_ids_hip

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




