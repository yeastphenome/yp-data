#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 26344761
paper_name = 'dobzinski_gerst_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[44]:


original_data = pd.read_excel('raw_data/table_S1.xlsx', sheet_name='Table 1', skiprows=1)


# In[45]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[46]:


original_data.head()


# In[47]:


original_data['orf'] = original_data['Systematic name'].astype(str)


# In[48]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[49]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[50]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[51]:


data_switch = {'+': 1, '±': 0.5, np.nan: 0}


# In[52]:


data_cols = ['Not targeted to vacuole under starvation conditions','Targeted to vacuole under nutrient-rich conditions']


# In[53]:


for c in data_cols:
    original_data[c] = original_data[c].apply(lambda x: data_switch[x])


# In[54]:


original_data.set_index('orf', inplace=True)


# In[55]:


original_data = original_data[data_cols].copy()


# In[56]:


original_data = original_data.groupby(original_data.index).mean()


# In[57]:


original_data.shape


# In[58]:


original_data.head()


# # Load & process tested strains

# In[59]:


tested = pd.read_excel('raw_data/KO_DAmP_ORFs.xlsx', sheet_name='Sheet1', skiprows=1)


# In[61]:


tested['orf'] = tested['ORF '].astype(str)


# In[62]:


tested['orf'] = clean_orf(tested['orf'])


# In[63]:


typo_fixes = {'YOLO57W':'YOL057W','YOLO62C':'YOL062C',
              'YJL206-A':'YJL206C-A','YLR287-A':'YLR287C-A','YBRF182C-A':'YBR182C-A'}


# In[64]:


for typo in typo_fixes.keys():
    tested.loc[tested['orf']==typo,'orf'] = typo_fixes[typo]


# In[65]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[66]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[67]:


tested = tested.loc[t,:]


# In[68]:


tested_orfs = np.unique(tested['orf'].values)


# In[69]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[70]:


# Removing the missing strains (they were tested as DAMP strains, not deletions?)
original_data.drop(index=missing, inplace=True)


# In[71]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[72]:


data = original_data.copy()


# In[73]:


dataset_ids = [16164, 16139]
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

data.head()


# # Normalize

# In[78]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[79]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[80]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

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




