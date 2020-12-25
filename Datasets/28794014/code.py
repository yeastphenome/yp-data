#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 28794014
paper_name = 'cohen_schuldiner_2017' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[46]:


original_data = pd.read_excel('raw_data/TableS1-S6.xlsx', sheet_name='S1', skiprows=1)


# In[47]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[48]:


original_data.head()


# In[49]:


original_data['orf'] = original_data['ORF'].astype(str)


# In[50]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[51]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[52]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[53]:


original_data['data'] = 0
original_data.loc[original_data['Grouping']=='Reduced clustering','data'] = -1
original_data.loc[original_data['Grouping']=='Enhanced clustering','data'] = 1


# In[54]:


original_data.set_index('orf', inplace=True)
original_data = original_data[['data']]


# In[55]:


original_data = original_data.groupby(original_data.index).mean()


# In[56]:


original_data.shape


# # Load & process tested strains

# In[57]:


tested = pd.read_excel('raw_data/KO_DAmP_ORFs.xlsx', sheet_name='Sheet1', skiprows=1)


# In[58]:


tested.head()


# In[59]:


tested.columns


# In[60]:


tested['orf'] = tested['ORF '].astype(str)
tested['orf'] = clean_orf(tested['orf'])


# In[61]:


typo_fix = {'YOLO57W':'YOL057W','YOLO62C':'YOL062C','YBRF182C-A':'YBR182C-A','YLR287-A':'YLR287C-A','YJL206-A':'YJL206C-A'}


# In[62]:


for t in typo_fix.keys():
    tested.loc[tested['orf']==t,'orf'] = typo_fix[t]


# In[63]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[64]:


t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[65]:


tested = tested.loc[t,:]


# In[66]:


tested_orfs = np.unique(tested['orf'].values)


# In[67]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[68]:


# Removing the missing strains (they were tested as DAMP strains, not deletions)
original_data.drop(index=missing, inplace=True)


# In[69]:


original_data.shape


# In[70]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# In[71]:


original_data.shape


# # Prepare the final dataset

# In[72]:


data = original_data[['data']].copy()


# In[73]:


dataset_ids = [15990]
datasets = datasets.reindex(index=dataset_ids)


# In[74]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[75]:


data.head()


# ## Subset to the genes currently in SGD

# In[76]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[77]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[78]:


data.head()


# # Normalize

# In[79]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[80]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[81]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[82]:


data_all.head()


# # Print out

# In[83]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[84]:


from IO.save_data_to_db3 import *


# In[85]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




