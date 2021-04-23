#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 32904421
paper_name = 'stenger_westermann_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# ## Data1

# In[5]:


original_data1 = pd.read_excel('raw_data/mic-07-234-s02.xlsx', sheet_name='Table S1')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))


# In[7]:


original_data1.head()


# In[9]:


original_data1['orf'] = original_data1['ORF'].astype(str)


# In[10]:


# Eliminate all white spaces & capitalize
original_data1['orf'] = clean_orf(original_data1['orf'])


# In[11]:


# Translate to ORFs 
original_data1['orf'] = translate_sc(original_data1['orf'], to='orf')


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orf'])
print(original_data1.loc[~t,])


# In[13]:


original_data1['data'] = original_data1['this study']


# In[14]:


data_switch = {'pet': -1, 'nd': np.nan, 0: 0}
original_data1['data'] = original_data1['data'].apply(lambda x: data_switch[x])


# In[15]:


original_data1['data'] = pd.to_numeric(original_data1['data'])


# In[16]:


original_data1.set_index('orf', inplace=True)


# In[17]:


original_data1 = original_data1[['data']].copy()


# In[18]:


original_data1 = original_data1.groupby(original_data1.index).mean()


# In[20]:


original_data1 = original_data1[original_data1['data'].notnull()]


# In[21]:


original_data1.shape


# ## Data2

# In[23]:


original_data2 = pd.read_csv('raw_data/Table2.txt', sep='\t', header=None)


# In[24]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[25]:


original_data2.head()


# In[26]:


original_data2['orf'] = original_data2[1].astype(str)


# In[27]:


# Eliminate all white spaces & capitalize
original_data2['orf'] = clean_orf(original_data2['orf'])


# In[28]:


# Translate to ORFs 
original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')


# In[29]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[30]:


original_data2['data'] = -1


# In[31]:


original_data2.set_index('orf', inplace=True)


# In[32]:


original_data2 = original_data2[['data']].copy()


# In[33]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[34]:


original_data2.shape


# In[35]:


tested_orfs = original_data1.index.values


# In[36]:


missing = [orf for orf in original_data2.index.values if orf not in tested_orfs]
missing


# In[37]:


original_data2 = original_data2.reindex(index=tested_orfs, fill_value=0)


# ## Data3

# In[50]:


original_data3 = pd.read_csv('raw_data/TableS3.txt', sep=' ')


# In[51]:


print('Original data dimensions: %d x %d' % (original_data3.shape))


# In[52]:


original_data3.head()


# In[53]:


original_data3['orf'] = original_data3['ORF'].astype(str)


# In[54]:


# Eliminate all white spaces & capitalize
original_data3['orf'] = clean_orf(original_data3['orf'])


# In[55]:


# Translate to ORFs 
original_data3['orf'] = translate_sc(original_data3['orf'], to='orf')


# In[56]:


# Make sure everything translated ok
t = looks_like_orf(original_data3['orf'])
print(original_data3.loc[~t,])


# In[57]:


original_data3['data'] = -1


# In[58]:


original_data3.set_index('orf', inplace=True)


# In[59]:


original_data3 = original_data3[['data']].copy()


# In[60]:


original_data3 = original_data3.groupby(original_data3.index).mean()


# In[61]:


original_data3.shape


# In[70]:


orfs1 = original_data1.index.values
orfs2 = original_data2.loc[original_data2['data']<0].index.values
tested_orfs = [orf for orf in orfs1 if orf not in orfs2]


# In[76]:


missing = [orf for orf in original_data3.index.values if orf not in tested_orfs]
missing


# In[77]:


original_data3 = original_data3.reindex(index=tested_orfs, fill_value=0)


# ## Data4

# In[92]:


original_data4 = pd.read_excel('raw_data/TableS5.xlsx', sheet_name='Table 1', skiprows=2)


# In[93]:


print('Original data dimensions: %d x %d' % (original_data4.shape))


# In[94]:


original_data4.head()


# In[95]:


original_data4['orf'] = original_data4['ORF'].astype(str)


# In[96]:


# Eliminate all white spaces & capitalize
original_data4['orf'] = clean_orf(original_data4['orf'])


# In[97]:


# Translate to ORFs 
original_data4['orf'] = translate_sc(original_data4['orf'], to='orf')


# In[98]:


# Make sure everything translated ok
t = looks_like_orf(original_data4['orf'])
print(original_data4.loc[~t,])


# In[99]:


original_data4 = original_data4.loc[t,:]


# In[100]:


original_data4['data'] = -1


# In[101]:


original_data4.set_index('orf', inplace=True)


# In[102]:


original_data4 = original_data4[['data']].copy()


# In[103]:


original_data4 = original_data4.groupby(original_data4.index).mean()


# In[104]:


original_data4.shape


# In[105]:


tested = pd.read_csv('raw_data/arg8_strains.txt', header=None)


# In[107]:


tested['orf'] = tested[0].astype(str)


# In[108]:


tested['orf'] = clean_orf(tested['orf'])


# In[109]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[110]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[111]:


tested_orfs = tested['orf'].unique()


# In[112]:


missing = [orf for orf in original_data4.index.values if orf not in tested_orfs]
missing


# In[113]:


original_data4 = original_data4.reindex(index=tested_orfs, fill_value=0)


# ## Pull all data together

# In[114]:


original_data = pd.concat([original_data1, original_data2, original_data3, original_data4], axis=1)


# In[117]:


original_data.index.name='orf'


# # Prepare the final dataset

# In[120]:


data = original_data.copy()


# In[121]:


dataset_ids = [21874, 21875, 21876, 21877]
datasets = datasets.reindex(index=dataset_ids)


# In[122]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[123]:


data.head()


# ## Subset to the genes currently in SGD

# In[124]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[125]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[126]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[127]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[128]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[129]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[130]:


from IO.save_data_to_db3 import *


# In[131]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




