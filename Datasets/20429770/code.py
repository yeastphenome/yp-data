#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 20429770
paper_name = 'mir_rashed_smith_2010' 


# In[49]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[50]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[24]:


original_data = pd.read_excel('raw_data/Ech sup Table work 5.xls', sheet_name='Sheet1', skiprows=5)


# In[25]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[26]:


original_data.head()


# In[27]:


original_data.columns.values


# In[29]:


cols_orfs = [ic for ic, c in enumerate(original_data.columns.values) if 'Gene ID' in c]
cols_orfs


# In[30]:


original_data_list = []
for c in cols_orfs:
    t = original_data.iloc[:,[c,c+1]].copy()
    t.columns = [0,1]
    t['orf'] = t[0].astype(str)
    t['orf'] = clean_orf(t['orf'])
    t['orf'] = translate_sc(t['orf'], to='orf')
    tx = looks_like_orf(t['orf'])
    print(t.loc[~tx,])
    t = t.loc[tx,:]
    t['data'] = -t[1].astype(float)
    t.set_index('orf', inplace=True)
    t = t[['data']].copy()
    t = t.groupby(t.index).mean()
    
    original_data_list.append(t)


# In[31]:


original_data = pd.concat(original_data_list, axis=1)


# In[32]:


original_data.index.name='orf'
original_data[original_data.isnull()] = 0


# In[45]:


col_names = ['SGEPH.Dark','SGEPH.UV','SGEPR.Dark','SGEPR.UV','SG1.Dark','SG1.UV',
             'SG7a.Dark','SG7a.UV','SGEPF.Dark','SGEPF.UV','SGEPLS.Dark',
             'SGEPLS.UV','SGEARa.Dark','SGEARa.UV','SGEARb.Dark','SGEARb.UV']


# In[47]:


original_data.columns = col_names


# In[48]:


original_data.head()


# # Prepare the final dataset

# In[51]:


data = original_data.copy()


# In[53]:


dataset_names = dict()
dataset_names['SG1.Dark'] = 16457
dataset_names['SG1.UV'] = 16575
dataset_names['SGEPR.Dark'] = 16576
dataset_names['SGEPR.UV'] = 16577
dataset_names['SGEPF.Dark'] = 16578
dataset_names['SGEPF.UV'] = 16579
dataset_names['SGEPH.Dark'] = 16580
dataset_names['SGEPH.UV'] = 16581
dataset_names['SGEPLS.Dark'] = 16582
dataset_names['SGEPLS.UV'] = 16583
dataset_names['SG7a.Dark'] = 16584
dataset_names['SG7a.UV'] = 16585
dataset_names['SGEARa.Dark'] = 16586
dataset_names['SGEARa.UV'] = 16587
dataset_names['SGEARb.Dark'] = 16588
dataset_names['SGEARb.UV'] = 16589


# In[54]:


dataset_ids = [dataset_names[c] for c in data.columns.values]
datasets = datasets.reindex(index=dataset_ids)


# In[55]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[56]:


data.head()


# ## Subset to the genes currently in SGD

# In[57]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[58]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[59]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[60]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[61]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[62]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[63]:


from IO.save_data_to_db3 import *


# In[64]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




