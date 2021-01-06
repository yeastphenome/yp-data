#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 23665193
paper_name = 'hirasawa_shimizu_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[95]:


original_data = pd.read_excel('raw_data/data_from_S1_PDF.xlsx', sheet_name='Mutants')


# In[96]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[97]:


original_data.head()


# In[98]:


t1 = pd.melt(original_data.iloc[:,[0,3,6,9]])
t2 = pd.melt(original_data.iloc[:,[1,4,7,10]])
t3 = pd.melt(original_data.iloc[:,[2,5,8,11]])


# In[99]:


original_data = pd.concat([t1,t2,t3], axis=1)


# In[100]:


original_data = original_data.loc[original_data.iloc[:,1].notnull(),:]


# In[101]:


original_data = original_data.iloc[:,[1,3,5]]
original_data.columns = ['genes','data1','data2']


# In[102]:


original_data['genes'] = original_data['genes'].astype(str)


# In[103]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[104]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['genes'].values, to='orf')


# In[105]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[106]:


original_data = original_data.loc[t,:]


# In[107]:


original_data.set_index('orf', inplace=True)


# In[108]:


original_data['data1'] = original_data['data1'].astype(float)
original_data['data2'] = original_data['data2'].astype(float)


# In[111]:


original_data['data'] = original_data[['data1','data2']].mean(axis=1)


# In[113]:


original_data = original_data[['data']].copy()


# In[114]:


original_data = original_data.groupby(original_data.index).mean()


# In[115]:


original_data.shape


# # Normalize by controls

# In[116]:


ctrl_data = pd.read_excel('raw_data/data_from_S1_PDF.xlsx', sheet_name='CTRL')


# In[117]:


print('Data dimensions: %d x %d' % (ctrl_data.shape))


# In[118]:


t2 = pd.melt(ctrl_data.iloc[:,[1,2,4,5,7,8,10,11]])


# In[119]:


ctrl_data = t2['value'].mean()


# In[120]:


original_data = original_data / ctrl_data


# # Prepare the final dataset

# In[122]:


data = original_data.copy()


# In[123]:


dataset_ids = [94]
datasets = datasets.reindex(index=dataset_ids)


# In[124]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[125]:


data.head()


# ## Subset to the genes currently in SGD

# In[126]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[127]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[128]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[129]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[130]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[131]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[132]:


from IO.save_data_to_db3 import *


# In[133]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




