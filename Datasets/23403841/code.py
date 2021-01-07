#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 23403841
paper_name = 'oconnor_vulpe_2012' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[7]:


original_data = {}
time = {}
for dt in np.arange(4)+1:
    original_data[dt] = pd.read_excel('raw_data/data sheet ' + str(dt) + '.xlsx', sheet_name='Sheet1', skiprows=1)
    time[dt] = original_data[dt].loc[0,'Treatment']
    original_data[dt].columns = ['genename','orf','82.5 mM','165 mM','330 mM']
    original_data[dt] = original_data[dt].drop(index=[0,1,2])
    
    original_data[dt]['orf'] = original_data[dt]['orf'].astype(str)
    original_data[dt]['orf'] = clean_orf(original_data[dt]['orf'])
    original_data[dt]['orf'] = translate_sc(original_data[dt]['orf'], to='orf')
    
    t = looks_like_orf(original_data[dt]['orf'])
    print(original_data[dt].loc[~t,])
    
    original_data[dt].set_index('orf', inplace=True)
    original_data[dt] = original_data[dt][['82.5 mM','165 mM','330 mM']].astype(float)
    original_data[dt] = original_data[dt].groupby(original_data[dt].index).mean()


# In[8]:


for dt in np.arange(4)+1:
    t = original_data[dt].copy()
    t.columns = [c+'_'+time[dt] for c in t.columns]
    if dt == 1:
        data_5g = t.copy()
    elif dt == 2:
        data_15g = t.copy()
    elif dt == 3:
        data_5g = pd.concat((data_5g, t), axis=0)
    elif dt == 4:
        data_15g = pd.concat((data_15g, t), axis=0)


# In[9]:


data_5g[data_5g.isnull()] = 0
data_15g[data_15g.isnull()] = 0


# In[10]:


data_5g = data_5g.astype(float)
data_15g = data_15g.astype(float)


# In[11]:


data_5g = data_5g.groupby(data_5g.index).mean()
data_15g = data_15g.groupby(data_15g.index).mean()


# In[12]:


data = data_5g.join(data_15g, how='outer')
data[data.isnull()] = 0


# In[13]:


data.shape


# In[14]:


data.head()


# # Prepare the final dataset

# In[15]:


dataset_ids = [16531,16529,16526,16530,16528,16527]
datasets = datasets.reindex(index=dataset_ids)


# In[16]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[17]:


data.head()


# ## Subset to the genes currently in SGD

# In[18]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[19]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[20]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[21]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[22]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[23]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[24]:


from IO.save_data_to_db3 import *


# In[25]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




