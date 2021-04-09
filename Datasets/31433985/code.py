#!/usr/bin/env python
# coding: utf-8

# In[3]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[4]:


paper_pmid = 31433985
paper_name = 'babazadeh_nystrom_2019' 


# In[5]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[6]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data (screen 1)

# In[7]:


original_data1 = pd.read_excel('raw_data/table_s1.xlsx', sheet_name='Table 1', skiprows=1)


# In[9]:


print('Original data dimensions: %d x %d' % (original_data1.shape))


# In[10]:


original_data1.head()


# In[16]:


original_data1['orf'] = original_data1['ORF'].astype(str).apply(lambda x: x.split(' ')[2] if len(x) > 3 else x)


# In[20]:


original_data1.loc[original_data1['orf']=='31','orf'] = 'VAC17'


# In[21]:


# Eliminate all white spaces & capitalize
original_data1['orf'] = clean_genename(original_data1['orf'])
original_data1['orf'] = clean_orf(original_data1['orf'])


# In[22]:


# Translate to ORFs 
original_data1['orf'] = translate_sc(original_data1['orf'], to='orf')


# In[23]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orf'])
print(original_data1.loc[~t,])


# In[24]:


original_data1['data'] = -1


# In[25]:


original_data1.set_index('orf', inplace=True)


# In[26]:


original_data1 = original_data1[['data']].copy()


# In[27]:


original_data1 = original_data1.groupby(original_data1.index).mean()


# In[28]:


original_data1.shape


# # Load & process the data (screen 2)

# In[49]:


original_data2 = pd.read_excel('raw_data/table_s2.xlsx', sheet_name='Sheet1')


# In[50]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[51]:


original_data2.head()


# In[52]:


original_data2['orf'] = original_data2['Systematic Name'].astype(str).apply(lambda x: x.split(' ')[2] if len(x) > 3 else x)


# In[53]:


# Eliminate all white spaces & capitalize
original_data2['orf'] = clean_orf(original_data2['orf'])


# In[54]:


# Translate to ORFs 
original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')


# In[55]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[56]:


# Remove essential data
original_data2['Gene/TS Allele'] = original_data2['Gene/TS Allele'].astype(str)
original_data2 = original_data2.loc[~original_data2['Gene/TS Allele'].str.contains('-')].copy()


# In[57]:


original_data2['data'] = pd.to_numeric(original_data2['Average Score'], errors='coerce')


# In[58]:


original_data2.set_index('orf', inplace=True)


# In[59]:


original_data2 = original_data2[['data']].copy()


# In[60]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[61]:


original_data2.shape


# # Merge

# In[70]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[73]:


original_data[original_data.isnull()] = 0


# In[74]:


original_data.head()


# # Load & process tested strains

# In[63]:


tested = pd.read_csv('raw_data/FG_array_genes.txt', sep='\t', header=None)


# In[64]:


tested.head()


# In[65]:


tested['orf'] = tested[0].astype(str)


# In[66]:


tested['orf'] = clean_orf(tested['orf'])


# In[67]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[68]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[69]:


tested_orfs = tested['orf'].unique()


# In[75]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[76]:


tested_orfs = list(tested_orfs) + missing


# In[77]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[78]:


data = original_data.copy()


# In[79]:


dataset_ids = [16664, 16663]
datasets = datasets.reindex(index=dataset_ids)


# In[80]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[81]:


data.head()


# ## Subset to the genes currently in SGD

# In[82]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[83]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[84]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[85]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[86]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[87]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[88]:


from IO.save_data_to_db3 import *


# In[90]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




