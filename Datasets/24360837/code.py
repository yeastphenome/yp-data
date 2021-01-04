#!/usr/bin/env python
# coding: utf-8

# In[80]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[81]:


paper_pmid = 24360837
paper_name = 'hoepfner_movva_2014' 


# In[82]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[83]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data - Benomyl

# In[84]:


original_data1 = pd.read_csv('large_files/raw_data/HOP_scores-benomyl.txt', sep='\t')
original_data2 = pd.read_csv('large_files/raw_data/HIP_scores-benomyl.txt', sep='\t')


# In[85]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[86]:


# Keep the sensitivity scores, not z-scores (z-score normalize each strain to its phenotype to all other compounds in the dataset)


# In[87]:


cols1 = [c for c in original_data1.columns.values if 'z-score' not in c]
cols2 = [c for c in original_data2.columns.values if 'z-score' not in c]


# In[88]:


original_data1 = original_data1.loc[:, cols1]
original_data2 = original_data2.loc[:, cols2]


# In[89]:


orf_col = 'Systematic Name'


# In[90]:


original_data1[orf_col] = original_data1[orf_col].astype(str)
original_data2[orf_col] = original_data2[orf_col].astype(str)


# In[91]:


# Eliminate all white spaces & capitalize
original_data1[orf_col] = clean_orf(original_data1[orf_col])
original_data2[orf_col] = clean_orf(original_data2[orf_col])


# In[92]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1[orf_col], to='orf')
original_data2['orfs'] = translate_sc(original_data2[orf_col], to='orf')


# In[93]:


original_data1.loc[original_data1['orfs'] == 'YBR160WAS','orfs'] = 'YBR160W'
original_data2.loc[original_data2['orfs'] == 'YBR160WAS','orfs'] = 'YBR160W'


# In[94]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[95]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[96]:


original_data1 = original_data1.loc[t,:]
original_data2 = original_data2.loc[t,:]


# In[97]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)
original_data1.index.name='orf'
original_data2.index.name='orf'


# In[98]:


original_data1['data'] = original_data1.mean(axis=1)
original_data2['data'] = original_data2.mean(axis=1)


# In[99]:


original_data = original_data1[['data']].join(original_data2[['data']], how='outer', lsuffix='_hop', rsuffix='_hip')


# In[100]:


original_data = original_data.groupby(original_data.index).mean()


# In[101]:


dataset_ids = [1087, 16622]
data_benomyl = original_data[['data_hop','data_hip']].copy()


# In[102]:


data_benomyl.columns = dataset_ids


# In[103]:


data_benomyl.head()


# # Load and process data -- all others

# In[104]:


original_data1 = pd.read_csv('large_files/raw_data/HOP_scores.txt', sep='\t')
original_data2 = pd.read_csv('large_files/raw_data/HIP_scores.txt', sep='\t')


# In[105]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[106]:


cols1 = [c for c in original_data1.columns.values if 'z-score' not in c]
cols2 = [c for c in original_data2.columns.values if 'z-score' not in c]


# In[107]:


original_data1 = original_data1.loc[:, cols1]
original_data2 = original_data2.loc[:, cols2]


# In[108]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


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
print(original_data1.loc[~t,])


# In[115]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[116]:


original_data1 = original_data1.loc[t,:]
original_data2 = original_data2.loc[t,:]


# In[117]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)
original_data1.index.name='orf'
original_data2.index.name='orf'


# In[118]:


original_data1.drop(columns=['Systematic Name'], inplace=True)


# In[119]:


original_data2.drop(columns=['Systematic Name'], inplace=True)


# ### Map dataset IDs to data columns

# In[120]:


compound_map = pd.read_csv('extras/type_cmb_dose_dataset.txt', sep='\t')


# In[121]:


compound_map.loc[compound_map.loc[:,'Dataset HOP']==1226]


# In[122]:


dt_ids = []
for s in original_data1.columns.values:
    s_parts = re.split(' |_',s)
    cmb = int(s_parts[4])
    dose = float(s_parts[5])
    
    dt = compound_map.loc[(compound_map['CMB'] == cmb) & (round(compound_map['Dose'],4) == round(dose,4))]
    if dt.shape[0] > 0:
        dataset_id = dt['Dataset HOP'].values[0]
    else:
        dataset_id = np.nan
    
    dt_ids.append(dataset_id)


# In[123]:


t = original_data1.drop(columns=original_data1.columns[np.isnan(np.array(dt_ids))])


# In[124]:


dt_ids = np.array(dt_ids)[~np.isnan(np.array(dt_ids))]


# In[125]:


dt_ids = dt_ids.astype(int)


# In[126]:


t.columns = dt_ids


# In[127]:


# Average values for duplicated (replicated) datasets
t = t.T
t = t.groupby(t.index).mean().T


# In[128]:


t.shape


# In[129]:


original_data1 = t.copy()


# In[130]:


dt_ids = []
for s in original_data2.columns.values:
    s_parts = re.split(' |_',s)
    cmb = int(s_parts[4])
    dose = float(s_parts[5])
    
    dt = compound_map.loc[(compound_map['CMB'] == cmb) & (round(compound_map['Dose'],4) == round(dose,4))]
    if dt.shape[0] > 0:
        dataset_id = dt['Dataset HIP'].values[0]
    else:
        dataset_id = np.nan
    
    dt_ids.append(dataset_id)


# In[131]:


t = original_data2.drop(columns=original_data2.columns[np.isnan(np.array(dt_ids))])


# In[132]:


dt_ids = np.array(dt_ids)[~np.isnan(np.array(dt_ids))]


# In[133]:


dt_ids = dt_ids.astype(int)


# In[134]:


t.columns = dt_ids


# In[135]:


# Average values for duplicated (replicated) datasets
t = t.T
t = t.groupby(t.index).mean().T


# In[136]:


t.shape


# In[137]:


original_data2 = t.copy()


# In[138]:


original_data2.shape


# ### Average and merge

# In[139]:


original_data1 = original_data1.groupby(original_data1.index).mean()


# In[140]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[141]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_hop', rsuffix='_hip')


# In[146]:


data_final = data_benomyl.join(original_data, how='outer', lsuffix='_benomyl', rsuffix='_other')


# In[147]:


data_final.shape


# In[148]:


data_benomyl.shape


# In[149]:


original_data.shape


# In[156]:


# Remove ORFs that are all NaNs
num_vals = data_final.notnull().sum(axis=1)


# In[157]:


data_final = data_final.loc[num_vals>0,:]


# In[158]:


data_final.shape


# # Prepare final dataset

# In[159]:


data = data_final.copy()


# In[160]:


dataset_ids = data_final.columns.values
datasets = datasets.reindex(index=dataset_ids)


# In[161]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[162]:


data.head()


# ## Subset to the genes currently in SGD

# In[163]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[164]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[165]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[166]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[167]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[168]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[169]:


from IO.save_data_to_db3 import *


# In[170]:


save_data_to_db(data_all, paper_pmid, delete=False)


# In[ ]:




