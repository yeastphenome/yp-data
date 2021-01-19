#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 18420932
paper_name = 'hillenmeyer_giaever_2008' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[43]:


original_data = pd.read_csv('raw_data/hom.ratio_result_nm.pub', sep='\t')


# In[44]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[45]:


original_data.head()


# In[46]:


original_data['orf'] = original_data['Orf'].apply(lambda x: x.split(':')[0])


# In[47]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[48]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[49]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[50]:


original_data.set_index('orf', inplace=True)


# In[51]:


original_data.drop(columns=['Orf'], inplace=True)


# In[52]:


original_data = original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[53]:


original_data = original_data.groupby(original_data.index).mean()


# In[54]:


original_data.shape


# ## Load dataset ids

# In[55]:


dt = pd.read_csv('extras/datasets.txt', sep='\t', header=None)


# In[56]:


dt[0] = pd.to_numeric(dt[0], errors='coerce')
dt.head()


# In[60]:


dt_ids = np.array([dt.loc[dt[1]==c,0].values[0] for c in original_data.columns.values])


# In[61]:


dt_ids[np.isnan(dt_ids)]


# In[62]:


original_data.columns = dt_ids


# In[66]:


original_data = original_data.loc[:,~np.isnan(original_data.columns)]


# In[68]:


original_data.columns = original_data.columns.values.astype(int)


# In[69]:


original_data.shape


# In[97]:


original_data = original_data.T
original_data = original_data.groupby(original_data.index).mean()
original_data = original_data.T


# In[98]:


original_data.shape


# # Load Het data

# In[71]:


original_data2 = pd.read_csv('raw_data/het.ratio_result_nm.pub', sep='\t')


# In[72]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[73]:


original_data2.head()


# In[74]:


original_data2['orf'] = original_data2['Orf'].apply(lambda x: x.split(':')[0])


# In[75]:


# Eliminate all white spaces & capitalize
original_data2['orf'] = clean_orf(original_data2['orf'])


# In[76]:


# Translate to ORFs 
original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')


# In[77]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[78]:


original_data2.set_index('orf', inplace=True)


# In[79]:


original_data2.drop(columns=['Orf'], inplace=True)


# In[80]:


original_data2 = original_data2.apply(pd.to_numeric, axis=1, errors='coerce')


# In[81]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[82]:


original_data2.shape


# ## Load dataset ids

# In[87]:


dt2 = pd.read_csv('extras/datasets_het.txt', sep='\t')


# In[88]:


dt2.shape


# In[89]:


dt2.head()


# In[90]:


dt2['Dataset'] = pd.to_numeric(dt2['Dataset'], errors='coerce')
dt2.head()


# In[92]:


dt_ids2 = np.array([dt2.loc[dt2['Name']==c,'Dataset'].values[0] for c in original_data2.columns.values])


# In[94]:


original_data2.columns = dt_ids2


# In[95]:


original_data2 = original_data2.loc[:,original_data2.columns>0]


# In[96]:


original_data2.shape


# In[99]:


original_data2 = original_data2.T
original_data2 = original_data2.groupby(original_data2.index).mean()
original_data2 = original_data2.T


# In[100]:


original_data2.shape


# # Merge

# In[101]:


original_data = original_data.join(original_data2, how='outer')


# In[102]:


original_data.shape


# In[103]:


original_data.head()


# In[117]:


# Taking the opposite because the original values are log2(control/treatment)
original_data = -original_data


# # Prepare the final dataset

# In[118]:


data = original_data.copy()


# In[119]:


dataset_ids = original_data.columns.values
datasets = datasets.reindex(index=dataset_ids)


# In[120]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[121]:


data.head()


# ## Subset to the genes currently in SGD

# In[122]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[123]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[124]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[125]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[126]:


vals = data.values
vals_norm = data_norm.values

vals_norm[np.isnan(vals)] = np.nan

data_norm = pd.DataFrame(index=data_norm.index, columns=data_norm.columns, data=vals_norm)


# In[127]:


data_all = data.join(data_norm)
data_all.head()


# # Print out

# In[128]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[129]:


# from IO.save_data_to_db3 import *


# In[130]:


# save_data_to_db(data_all, paper_pmid)

