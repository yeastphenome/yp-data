#!/usr/bin/env python
# coding: utf-8

# In[32]:


get_ipython().run_line_magic('run', '../yp_utils.py')

import matplotlib.pyplot as plt


# # Initial setup

# In[2]:


paper_pmid = 32994210
paper_name = 'stjohn_fasullo_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[18]:


files = ['GSE129699_C8B0CANXX_Results__No_CYP_uptag_only__DMSO_AFB1.xls', 
         'GSE129699_C8B0CANXX_Results__CYP1A2_uptag_only__DMSO_AFB1.xls', 
         'GSE129699_C8B0CANXX_Results__CYP1A2nat214_uptag_only__DMSO_AFB1.xls']


# In[34]:


original_data = []
for f in files:
    path_to_file = 'raw_data/' + f
    data = pd.read_excel(path_to_file, sheet_name=None)
    
    sheet_name = list(data.keys())[0]
    original_data.append(data[sheet_name])


# In[38]:


original_data2 = []
for df in original_data:
    df['orfs'] = df['gene_id'].astype(str)
    df['orfs'] = clean_orf(df['orfs'])
    df['orfs'] = translate_sc(df['orfs'], to='orf')
    t = looks_like_orf(df['orfs'])
    df = df.loc[t,:]
    df = df.groupby('orfs').mean()
    df.index.name='orf'
    df['data'] = df['m.value']
    
    # For some reason, essential genes are present in the list, but they all have low values
    # So, it's probably an effect of mapping dictionary. Will remove
    ess = is_essential(df.index.values)
    df.drop(index=df.index.values[ess.values], inplace=True)
    
    print(df.shape)
    original_data2.append(df)


# In[41]:


data = pd.concat([df['data'] for df in original_data2], axis=1)


# In[42]:


data.shape


# In[43]:


data.head()


# # Prepare the final dataset

# In[44]:


dataset_ids = [16686,16647,16687]
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


# In[49]:


data.head()


# # Normalize

# In[50]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[51]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[52]:


data_norm[data.isnull()] = np.nan


# In[53]:


data_all = data.join(data_norm)


# In[54]:


data_all.head()


# In[58]:


data_all.sort_values(by=(16686,'value'))


# # Print out

# In[55]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[59]:


from IO.save_data_to_db3 import *


# In[60]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




