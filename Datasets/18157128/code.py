#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 18157128
paper_name = 'delneri_oliver_2008' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[91]:


original_data = pd.read_excel('raw_data/Table 4s_NewForPublishing.xlsx', sheet_name='CL, NL, PL & GJ  data', skiprows=10)


# In[92]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[93]:


original_data.head()


# In[94]:


orf_cols = [c for c in original_data.columns if 'ORF' in c]
data_cols = [c for c in original_data.columns if 'Growth Rate' in c]


# In[95]:


original_data_list = []
for o,d in zip(orf_cols, data_cols):
    t = original_data.loc[:,[o,d]].copy()
    t.columns = ['orf','data']
    original_data_list.append(t)


# In[96]:


original_data_list2 = []
for d in np.arange(4):
    df = pd.concat(original_data_list[d*2:d*2+1], axis=0)
    
    df['orf'] = df['orf'].astype(str)
    df['orf'] = clean_orf(df['orf'])
    df['orf'] = translate_sc(df['orf'].values, to='orf')
    
    # Make sure everything translated ok
    t = looks_like_orf(df['orf'])
    print(df.loc[~t,])
    df = df.loc[t,:]
    
    df.set_index('orf', inplace=True)
    df = df[['data']].copy()
    df['data'] = pd.to_numeric(df['data'], errors='coerce')
    
    df = df.groupby(df.index).mean()
    print(df.shape)
    
    original_data_list2.append(df)


# In[97]:


original_data = pd.concat(original_data_list2, axis=1)


# In[98]:


original_data.head()


# In[99]:


original_data.index.name = 'orf'


# In[100]:


original_data.shape


# # Prepare the final dataset

# In[101]:


data = original_data.copy()


# In[102]:


dataset_ids = [11813,11815,11816,11814]
datasets = datasets.reindex(index=dataset_ids)


# In[103]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[104]:


data.head()


# ## Subset to the genes currently in SGD

# In[105]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[106]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[107]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[108]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[109]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[110]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[111]:


from IO.save_data_to_db3 import *


# In[112]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




