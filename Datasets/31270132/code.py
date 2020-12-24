#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 31270132
paper_name = 'hoffert_strome_2019' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[46]:


original_data = pd.read_excel('raw_data/Table_S1.xlsx', sheet_name='Sheet1', skiprows=2)


# In[47]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[48]:


original_data['ORF name'] = original_data['ORF name'].astype(str)


# In[49]:


# Eliminate all white spaces & capitalize
original_data['ORF name'] = clean_orf(original_data['ORF name'])


# In[50]:


# Translate to ORFs 
original_data['ORF name'] = translate_sc(original_data['ORF name'], to='orf')


# In[51]:


typo_fixes = {'YCLO51W': 'YCL051W','YHR139C-': 'YHR139C-A','YGR122C-': 'YGR122C-A'}


# In[52]:


for orf in typo_fixes.keys():
    original_data.loc[original_data['ORF name']==orf,'ORF name'] = typo_fixes[orf]


# In[53]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF name'])


# In[54]:


# Remove the 1's at the end of certain ORFs
for orf in original_data.loc[~t,'ORF name'].values:
    new_orf = orf.rstrip('1')
    original_data.loc[original_data['ORF name']==orf,'ORF name'] = new_orf


# In[55]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF name'])
print(original_data.loc[~t,])


# In[56]:


original_data = original_data.loc[t,]


# In[57]:


original_data.head()


# In[58]:


data_replacements = {0: 0, '+': 1, '++': 2, '+++': 3, '++++': 4}


# In[59]:


original_data = original_data[['ORF name','A','B','A.1','B.1','A.2','B.2','A.3','B.3']]


# In[60]:


original_data.set_index('ORF name', inplace=True)
original_data.index.name='orf'


# In[61]:


for c in original_data.columns:
    original_data[c+'_num'] = original_data[c].apply(lambda x: data_replacements[x] if x in data_replacements else np.nan)


# In[62]:


original_data['data'] = original_data[['A_num','B_num','A.1_num','B.1_num','A.2_num','B.2_num','A.3_num','B.3_num']].sum(axis=1)


# In[63]:


original_data.shape


# In[64]:


original_data['num_vals'] = original_data[['A_num','B_num','A.1_num','B.1_num','A.2_num','B.2_num','A.3_num','B.3_num']].apply(lambda x: ~np.isnan(x)).sum(axis=1)


# In[65]:


original_data = original_data.loc[original_data['num_vals']>0,]


# In[66]:


original_data = original_data[['data']].copy()


# In[67]:


original_data = original_data.groupby(original_data.index).mean()


# In[68]:


original_data.shape


# # Prepare the final dataset

# In[69]:


data = original_data[['data']].copy()


# In[70]:


dataset_ids = [16548]
datasets = datasets.reindex(index=dataset_ids)


# In[71]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[72]:


data.head()


# ## Subset to the genes currently in SGD

# In[73]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[74]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[75]:


data.head()


# # Normalize

# In[76]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[77]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[78]:


data_norm[data.isnull()] = np.nan


# In[79]:


data_all = data.join(data_norm)


# In[80]:


data_all.head()


# # Print out

# In[81]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[82]:


from IO.save_data_to_db3 import *


# In[83]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




