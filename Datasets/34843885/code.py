#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 34843885
paper_name = 'guan_zhang_2022' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


sheets = ['4-NQO','FA','DCA']
original_data_list = []
for s in sheets:
    original_data_list.append(pd.read_excel('raw_data/1-s2.0-S0887233321002034-mmc1.xlsx', sheet_name=s, skiprows=1))


# In[6]:


print('Original data dimensions: %d x %d' % (original_data_list[0].shape))


# In[7]:


original_data_list[0].head()


# In[8]:


for s in np.arange(3):
    original_data_list[s]['orf'] = original_data_list[s]['Gene ORF'].astype(str)
    original_data_list[s]['orf'] = clean_orf(original_data_list[s]['orf'])
    original_data_list[s]['orf'] = translate_sc(original_data_list[s]['orf'], to='orf')
    t = looks_like_orf(original_data_list[s]['orf'])
    print(original_data_list[s].loc[~t,])
    original_data_list[s].set_index('orf', inplace=True)
    original_data_list[s].drop(columns=['Gene ORF'], inplace=True)
    for c in original_data_list[s].columns:
        original_data_list[s][c] = pd.to_numeric(original_data_list[s][c], errors='coerce')


# In[9]:


original_data_list[0].head()


# In[38]:


controls = ['DMSO','DMSO','H2O']
original_data_list2 = []
for s in np.arange(3):
    exps = ['_'.join(c.split('_')[:-1]) for c in original_data_list[s].columns]
    t = original_data_list[s].groupby(by=exps, axis=1).mean()
    
    # Remove mutants with low counts (less than 10) from any control group (DMSO/H2O)
    t = t.loc[t[controls[s]]>10,:]
    
    # Normalize by control
    t = t.div(t[controls[s]], axis=0)
    t.drop(columns=controls[s], inplace=True)
    t = t.groupby(t.index).mean()
    
    # Remove NaNs
    n = np.sum(~np.isnan(t.values) & ~np.isinf(t.values), axis=1)
    t = t.loc[n>0,:]
    t[np.isinf(t)] = np.nan
    
    original_data_list2.append(t)


# In[39]:


original_data = pd.concat(original_data_list2, axis=1)


# In[40]:


original_data.shape


# In[41]:


original_data.head()


# In[42]:


original_data.rename_axis(index='orf', inplace=True)


# # Prepare the final dataset

# In[43]:


data = original_data.copy()


# In[44]:


dataset_ids = [22101, 22100, 22099, 22098, 22097, 22105, 22104, 22103, 22109, 22108, 22107]
datasets = datasets.reindex(index=dataset_ids)


# In[45]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[46]:


data.head()


# ## Subset to the genes currently in SGD

# In[47]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[48]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[49]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[50]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[51]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[52]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[53]:


from IO.save_data_to_db3 import *


# In[54]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




