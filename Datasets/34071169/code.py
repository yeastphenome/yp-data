#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 34071169
paper_name = 'kipanga_luyten_2021' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[57]:


dt1 = pd.read_excel('raw_data/TablesS1-S2.xlsx', sheet_name='S1', header=None, index=0)
dt2 = pd.read_excel('raw_data/TablesS1-S2.xlsx', sheet_name='S2', header=None, index=0)


# In[61]:


dt1.set_index(0, inplace=True)
dt2.set_index(0, inplace=True)


# In[64]:


original_data = dt1.join(dt2, how='outer', lsuffix='_s1', rsuffix='_s2')


# In[65]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[66]:


original_data.head()


# In[67]:


original_data['genes'] = original_data.index.values.astype(str)


# In[68]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[73]:


original_data['genes'] = original_data['genes'].astype(str)


# In[75]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['genes'].values, to='orf')


# In[77]:


original_data.loc[original_data['genes']=='ATG42','orf'] = 'YBR139W'


# In[78]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[12]:





# In[79]:


original_data.set_index('orf', inplace=True)


# In[80]:


original_data = original_data[['1_s1','1_s2']].copy()


# In[81]:


original_data = original_data.groupby(original_data.index).mean()


# In[82]:


original_data.shape


# In[83]:


original_data.head()


# # Load & process tested strains

# In[86]:


tested = pd.read_excel('raw_data/screening library.xlsx', sheet_name='List Of Strains', skiprows=3, header=None)


# In[87]:


tested.head()


# In[88]:


tested['orf'] = tested[1].astype(str)


# In[89]:


tested['orf'] = clean_orf(tested['orf'])


# In[90]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[91]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[92]:


tested = tested.loc[t,:]


# In[93]:


tested_orfs = tested['orf'].unique()


# In[94]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[95]:


original_data = original_data.reindex(index=tested_orfs, fill_value=3.5)


# In[97]:


original_data[original_data.isnull()] = 3.5


# # Prepare the final dataset

# In[98]:


data = original_data.copy()


# In[99]:


dataset_ids = [21950,21951]
datasets = datasets.reindex(index=dataset_ids)


# In[100]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[101]:


data.head()


# ## Subset to the genes currently in SGD

# In[102]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[103]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[104]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[105]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[106]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[107]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[108]:


from IO.save_data_to_db3 import *


# In[109]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




