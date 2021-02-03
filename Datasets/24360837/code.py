#!/usr/bin/env python
# coding: utf-8

# In[219]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[220]:


paper_pmid = 24360837
paper_name = 'hoepfner_movva_2014' 


# In[221]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[222]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data - Benomyl

# In[104]:


original_data1 = pd.read_csv('large_files/raw_data/HOP_scores-benomyl.txt', sep='\t')
original_data2 = pd.read_csv('large_files/raw_data/HIP_scores-benomyl.txt', sep='\t')


# In[105]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[106]:


# Keep the sensitivity scores, not z-scores (z-score normalize each strain to its phenotype to all other compounds in the dataset)


# In[107]:


cols1 = [c for c in original_data1.columns.values if 'z-score' not in c]
cols2 = [c for c in original_data2.columns.values if 'z-score' not in c]


# In[108]:


original_data1 = original_data1.loc[:, cols1]
original_data2 = original_data2.loc[:, cols2]


# In[109]:


orf_col = 'Systematic Name'


# In[110]:


original_data1[orf_col] = original_data1[orf_col].astype(str)
original_data2[orf_col] = original_data2[orf_col].astype(str)


# In[111]:


# Eliminate all white spaces & capitalize
original_data1[orf_col] = clean_orf(original_data1[orf_col])
original_data2[orf_col] = clean_orf(original_data2[orf_col])


# In[112]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1[orf_col], to='orf')
original_data2['orfs'] = translate_sc(original_data2[orf_col], to='orf')


# In[113]:


original_data1.loc[original_data1['orfs'] == 'YBR160WAS','orfs'] = 'YBR160W'
original_data2.loc[original_data2['orfs'] == 'YBR160WAS','orfs'] = 'YBR160W'


# In[114]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
# print(original_data1.loc[~t,])


# In[115]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
# print(original_data2.loc[~t,])


# In[116]:


original_data1 = original_data1.loc[t,:]
original_data2 = original_data2.loc[t,:]


# In[117]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)
original_data1.index.name='orf'
original_data2.index.name='orf'


# In[118]:


original_data1['data'] = original_data1.mean(axis=1)
original_data2['data'] = original_data2.mean(axis=1)


# In[119]:


original_data = original_data1[['data']].join(original_data2[['data']], how='outer', lsuffix='_hop', rsuffix='_hip')


# In[120]:


original_data = original_data.groupby(original_data.index).mean()


# In[121]:


dataset_ids = [1087, 16622]
data_benomyl = original_data[['data_hop','data_hip']].copy()


# In[122]:


data_benomyl.columns = dataset_ids


# In[123]:


data_benomyl.head()


# # Load and process data -- all others

# In[232]:


original_data1 = pd.read_csv('large_files/raw_data/HOP_scores.txt', sep='\t')
original_data2 = pd.read_csv('large_files/raw_data/HIP_scores.txt', sep='\t')


# In[239]:


original_data1.set_index('Systematic Name', inplace=True)


# In[240]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[245]:


random_rows = np.random.choice(original_data1.index, 5)
random_cols = np.random.choice(original_data1.columns, 5)
original_data1.loc[random_rows, random_cols]


# In[227]:


cols1 = [c for c in original_data1.columns.values if 'z-score' not in c]
cols2 = [c for c in original_data2.columns.values if 'z-score' not in c]


# In[228]:


original_data1 = original_data1.loc[:, cols1]
original_data2 = original_data2.loc[:, cols2]


# In[229]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[129]:


orf_col = 'Systematic Name'


# In[130]:


original_data1[orf_col] = original_data1[orf_col].astype(str)
original_data2[orf_col] = original_data2[orf_col].astype(str)


# In[131]:


# Eliminate all white spaces & capitalize
original_data1[orf_col] = clean_orf(original_data1[orf_col])
original_data2[orf_col] = clean_orf(original_data2[orf_col])


# In[132]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1[orf_col], to='orf')
original_data2['orfs'] = translate_sc(original_data2[orf_col], to='orf')


# In[133]:


original_data1.loc[original_data1['orfs'] == 'YBR160WAS','orfs'] = 'YBR160W'
original_data2.loc[original_data2['orfs'] == 'YBR160WAS','orfs'] = 'YBR160W'


# In[134]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
# print(original_data1.loc[~t,])


# In[135]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
# print(original_data2.loc[~t,])


# In[136]:


original_data1 = original_data1.loc[t,:]
original_data2 = original_data2.loc[t,:]


# In[137]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)
original_data1.index.name='orf'
original_data2.index.name='orf'


# In[138]:


original_data1.drop(columns=['Systematic Name'], inplace=True)


# In[139]:


original_data2.drop(columns=['Systematic Name'], inplace=True)


# In[140]:


original_data1 = original_data1.groupby(original_data1.index).mean()


# In[141]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# ### Map data columns to dataset_ids

# In[142]:


dt = pd.read_csv('extras/datasets_name_to_id.txt', sep='\t')


# In[143]:


dt.head()


# In[144]:


dt.set_index('name', inplace=True)


# In[145]:


dt1 = dt.reindex(index=original_data1.columns.values)


# In[146]:


dt2 = dt.reindex(index=original_data2.columns.values)


# In[147]:


dt2.head()


# In[148]:


original_data1.columns = dt1['dataset'].values
original_data1 = original_data1.T
original_data1 = original_data1.groupby(original_data1.index).mean()
original_data1 = original_data1.T
original_data1.shape


# In[149]:


original_data2.columns = dt2['dataset'].values
original_data2 = original_data2.T
original_data2 = original_data2.groupby(original_data2.index).mean()
original_data2 = original_data2.T
original_data2.shape


# ### Merge

# In[150]:


original_data = original_data1.join(original_data2, how='outer')


# In[151]:


original_data.shape


# In[152]:


original_data_final = data_benomyl.join(original_data, how='outer', lsuffix='_benomyl', rsuffix='_other')


# In[153]:


original_data_final.shape


# In[154]:


data_benomyl.shape


# In[155]:


original_data.shape


# In[156]:


# Remove ORFs that are all NaNs
num_vals = original_data_final.notnull().sum(axis=1)


# In[157]:


original_data_final = original_data_final.loc[num_vals>0,:]


# In[158]:


original_data_final.shape


# In[159]:


original_data_final.head()


# # Prepare final dataset

# In[160]:


data = original_data_final.copy()


# In[161]:


dataset_ids = original_data_final.columns.values
datasets = datasets.reindex(index=dataset_ids)


# In[162]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[163]:


data.head()


# ## Subset to the genes currently in SGD

# In[164]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[165]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[166]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[167]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[168]:


data_vals = data.values
data_norm_vals = data_norm.values

data_norm_vals[np.isnan(data_vals)] = np.nan

data_norm = pd.DataFrame(index=data_norm.index, columns=data_norm.columns, data=data_norm_vals)


# In[169]:


data_all = data.join(data_norm)
data_all.head()


# # Print out

# In[71]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[170]:


# from IO.save_data_to_db3 import *


# In[171]:


# save_data_to_db(data_all, paper_pmid, delete=True)


# In[ ]:




