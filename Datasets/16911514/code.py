#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 16911514
paper_name = 'kawahata_iefuji_2006' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[42]:


sheets = ['Resistant','Sensitive']

original_data_list = []
for ixs, s in enumerate(sheets):
    original_data = pd.read_excel('raw_data/hits.xlsx', sheet_name=s)
    print('Original data dimensions: %d x %d' % (original_data.shape))
    print(original_data.head())
    
    original_data1 = original_data.iloc[:,[0,2,3,4]].copy()
    original_data1.columns = ['orf','lactic','acetic','hydrochloric']
    
    original_data2 = original_data.iloc[:,[5,7,8,9]].copy()
    original_data2.columns = ['orf','lactic','acetic','hydrochloric']
    
    original_data = pd.concat([original_data1, original_data2], axis=0)
    
    original_data['orf'] = original_data['orf'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'].values, to='orf')
    
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    
    original_data = original_data.loc[t,:]
    
    original_data.set_index('orf', inplace=True)
    
    for c in original_data.columns:
        original_data[c] = original_data[c].apply(lambda x: 0 if str(x) == 'nan' else np.power(-1, ixs) * len(x))
    
    original_data = original_data.groupby(original_data.index).mean()
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[43]:


original_data = pd.concat(original_data_list, axis=1)


# In[44]:


original_data.shape


# In[45]:


original_data = original_data.groupby(original_data.index).mean()


# In[46]:


original_data.shape


# In[47]:


original_data.index.name = 'orf'


# In[48]:


original_data[original_data.isnull()] = 0


# In[49]:


original_data.head()


# # Prepare the final dataset

# In[64]:


data = original_data.copy()


# In[65]:


dataset_ids = [422, 420, 418, 177, 421, 419]
datasets = datasets.reindex(index=dataset_ids)


# In[66]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[67]:


data.head()


# ## Subset to the genes currently in SGD

# In[68]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[69]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[70]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[71]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[72]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[73]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[74]:


from IO.save_data_to_db3 import *


# In[75]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




