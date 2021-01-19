#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 18245339
paper_name = 'abe_minegishi_2008' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[134]:


original_data = pd.read_excel('raw_data/Table1_Abe_Genetics.xlsx', sheet_name='Table 1 (2)', skiprows=5)


# In[135]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[136]:


original_data.head()


# In[137]:


original_data['orf'] = original_data['Unnamed: 2'].astype(str)


# In[138]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[139]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[140]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[141]:


original_data = original_data.loc[t,:]


# In[142]:


original_data.set_index('orf', inplace=True)


# In[143]:


# Data originally reported as percent relative to WT (100%). 
# So we're scaling back to fraction and shifting by 1, so that lower percentages correspond to the most negative values
original_data = original_data.iloc[:,[25,30]]/100 - 1


# In[144]:


original_data = original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[145]:


original_data = original_data.groupby(original_data.index).mean()


# In[146]:


original_data.shape


# # Load & process tested strains

# In[147]:


tested = pd.read_excel('raw_data/mat_alpha_041902.xlsx', sheet_name='mat_alpha_041902.txt', skiprows=1)


# In[148]:


tested.head()


# In[149]:


tested['orf'] = tested['ORF name'].astype(str)


# In[150]:


tested['orf'] = clean_orf(tested['orf'])


# In[151]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[152]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[153]:


tested = tested.loc[t,:]


# In[154]:


tested_orfs = tested['orf'].unique()


# In[155]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[156]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[157]:


data = original_data.copy()


# In[158]:


dataset_ids = [537,538]
datasets = datasets.reindex(index=dataset_ids)


# In[159]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[160]:


data.head()


# ## Subset to the genes currently in SGD

# In[161]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[162]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[163]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[164]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[165]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[166]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[167]:


from IO.save_data_to_db3 import *


# In[168]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




