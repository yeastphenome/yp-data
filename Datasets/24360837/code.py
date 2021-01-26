#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 24360837
paper_name = 'hoepfner_movva_2014' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data - Benomyl

# In[5]:


original_data1 = pd.read_csv('large_files/raw_data/HOP_scores-benomyl.txt', sep='\t')
original_data2 = pd.read_csv('large_files/raw_data/HIP_scores-benomyl.txt', sep='\t')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[7]:


# Keep the sensitivity scores, not z-scores (z-score normalize each strain to its phenotype to all other compounds in the dataset)


# In[8]:


cols1 = [c for c in original_data1.columns.values if 'z-score' not in c]
cols2 = [c for c in original_data2.columns.values if 'z-score' not in c]


# In[9]:


original_data1 = original_data1.loc[:, cols1]
original_data2 = original_data2.loc[:, cols2]


# In[10]:


orf_col = 'Systematic Name'


# In[11]:


original_data1[orf_col] = original_data1[orf_col].astype(str)
original_data2[orf_col] = original_data2[orf_col].astype(str)


# In[12]:


# Eliminate all white spaces & capitalize
original_data1[orf_col] = clean_orf(original_data1[orf_col])
original_data2[orf_col] = clean_orf(original_data2[orf_col])


# In[13]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1[orf_col], to='orf')
original_data2['orfs'] = translate_sc(original_data2[orf_col], to='orf')


# In[14]:


original_data1.loc[original_data1['orfs'] == 'YBR160WAS','orfs'] = 'YBR160W'
original_data2.loc[original_data2['orfs'] == 'YBR160WAS','orfs'] = 'YBR160W'


# In[15]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[16]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[17]:


original_data1 = original_data1.loc[t,:]
original_data2 = original_data2.loc[t,:]


# In[18]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)
original_data1.index.name='orf'
original_data2.index.name='orf'


# In[19]:


original_data1['data'] = original_data1.mean(axis=1)
original_data2['data'] = original_data2.mean(axis=1)


# In[20]:


original_data = original_data1[['data']].join(original_data2[['data']], how='outer', lsuffix='_hop', rsuffix='_hip')


# In[21]:


original_data = original_data.groupby(original_data.index).mean()


# In[22]:


dataset_ids = [1087, 16622]
data_benomyl = original_data[['data_hop','data_hip']].copy()


# In[23]:


data_benomyl.columns = dataset_ids


# In[24]:


data_benomyl.head()


# # Load and process data -- all others

# In[25]:


original_data1 = pd.read_csv('large_files/raw_data/HOP_scores.txt', sep='\t')
original_data2 = pd.read_csv('large_files/raw_data/HIP_scores.txt', sep='\t')


# In[26]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[27]:


cols1 = [c for c in original_data1.columns.values if 'z-score' not in c]
cols2 = [c for c in original_data2.columns.values if 'z-score' not in c]


# In[28]:


original_data1 = original_data1.loc[:, cols1]
original_data2 = original_data2.loc[:, cols2]


# In[29]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[30]:


orf_col = 'Systematic Name'


# In[31]:


original_data1[orf_col] = original_data1[orf_col].astype(str)
original_data2[orf_col] = original_data2[orf_col].astype(str)


# In[32]:


# Eliminate all white spaces & capitalize
original_data1[orf_col] = clean_orf(original_data1[orf_col])
original_data2[orf_col] = clean_orf(original_data2[orf_col])


# In[33]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1[orf_col], to='orf')
original_data2['orfs'] = translate_sc(original_data2[orf_col], to='orf')


# In[34]:


original_data1.loc[original_data1['orfs'] == 'YBR160WAS','orfs'] = 'YBR160W'
original_data2.loc[original_data2['orfs'] == 'YBR160WAS','orfs'] = 'YBR160W'


# In[35]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[36]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[37]:


original_data1 = original_data1.loc[t,:]
original_data2 = original_data2.loc[t,:]


# In[38]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)
original_data1.index.name='orf'
original_data2.index.name='orf'


# In[39]:


original_data1.drop(columns=['Systematic Name'], inplace=True)


# In[40]:


original_data2.drop(columns=['Systematic Name'], inplace=True)


# ### Map dataset IDs to data columns

# In[41]:


compound_map = pd.read_csv('extras/type_cmb_dose_dataset.txt', sep='\t')


# In[42]:


compound_map.loc[compound_map.loc[:,'Dataset HOP']==1226]


# In[43]:


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


# In[44]:


t = original_data1.drop(columns=original_data1.columns[np.isnan(np.array(dt_ids))])


# In[45]:


dt_ids = np.array(dt_ids)[~np.isnan(np.array(dt_ids))]


# In[46]:


dt_ids = dt_ids.astype(int)


# In[47]:


t.columns = dt_ids


# In[48]:


# Average values for duplicated (replicated) datasets
t = t.T
t = t.groupby(t.index).mean().T


# In[49]:


t.shape


# In[50]:


original_data1 = t.copy()


# In[51]:


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


# In[52]:


t = original_data2.drop(columns=original_data2.columns[np.isnan(np.array(dt_ids))])


# In[53]:


dt_ids = np.array(dt_ids)[~np.isnan(np.array(dt_ids))]


# In[54]:


dt_ids = dt_ids.astype(int)


# In[55]:


t.columns = dt_ids


# In[56]:


# Average values for duplicated (replicated) datasets
t = t.T
t = t.groupby(t.index).mean().T


# In[57]:


t.shape


# In[58]:


original_data2 = t.copy()


# In[59]:


original_data2.shape


# ### Average and merge

# In[60]:


original_data1 = original_data1.groupby(original_data1.index).mean()


# In[61]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[62]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_hop', rsuffix='_hip')


# In[63]:


data_final = data_benomyl.join(original_data, how='outer', lsuffix='_benomyl', rsuffix='_other')


# In[64]:


data_final.shape


# In[65]:


data_benomyl.shape


# In[66]:


original_data.shape


# In[67]:


# Remove ORFs that are all NaNs
num_vals = data_final.notnull().sum(axis=1)


# In[68]:


data_final = data_final.loc[num_vals>0,:]


# In[69]:


data_final.shape


# In[70]:


data_final.head()


# # Prepare final dataset

# In[71]:


data = data_final.copy()


# In[72]:


dataset_ids = data_final.columns.values
datasets = datasets.reindex(index=dataset_ids)


# In[73]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[74]:


data.head()


# ## Subset to the genes currently in SGD

# In[75]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[76]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[77]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[78]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[80]:


data_vals = data.values
data_norm_vals = data_norm.values

data_norm_vals[np.isnan(data_vals)] = np.nan

data_norm = pd.DataFrame(index=data_norm.index, columns=data_norm.columns, data=data_norm_vals)


# In[81]:


data_all = data.join(data_norm)
data_all.head()


# # Print out

# In[82]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[83]:


from IO.save_data_to_db3 import *


# In[84]:


save_data_to_db(data_all, paper_pmid, delete=True)


# In[ ]:




