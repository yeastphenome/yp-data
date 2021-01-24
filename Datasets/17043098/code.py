#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 17043098
paper_name = 'doostzadeh_langston_2007' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[47]:


original_data = pd.read_excel('raw_data/kfl131supp.xls', sheet_name='ToxSci Supplementary Data files', skiprows=1)


# In[48]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[49]:


original_data['strain'] = original_data['strain'].astype(str)


# In[50]:


mpp_columns = ['z_result_nq:04_10_28_19:mpp+:250:ug/ml::::20:hom_09_02',
              'z_result_nq:04_10_28_25:mpp+:250:ug/ml::::20:hom_09_02',
              'z_result_nq:04_11_04_06:mpp+:250:ug/ml::::20:hom_09_02']

paraquat_columns = ['z_result_nq:04_10_28_21:paraquat:5000:uM::::20:hom_09_02',
                   'z_result_nq:04_10_28_27:paraquat:5000:uM::::20:hom_09_02',
                   'z_result_nq:04_11_04_08:paraquat:5000:uM::::20:hom_09_02']


# In[51]:


# Extract ORF from string
orfs = original_data['strain'].apply(lambda x: x.split(':')[0])


# In[52]:


original_data['orfs'] = orfs


# In[53]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[54]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[55]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[56]:


original_data['mpp'] = original_data[mpp_columns].apply(pd.to_numeric, axis=1, errors='coerce').mean(axis=1)


# In[57]:


original_data['paraquat'] = original_data[paraquat_columns].apply(pd.to_numeric, axis=1, errors='coerce').mean(axis=1)


# In[58]:


original_data.set_index('orfs', inplace=True)


# In[59]:


original_data.index.name='orf'


# In[60]:


original_data = -original_data[['mpp','paraquat']].copy()


# In[61]:


original_data = original_data.groupby(original_data.index).mean()


# In[62]:


original_data.shape


# # Prepare the final dataset

# In[63]:


data = original_data.copy()


# In[64]:


dataset_ids = [16616,16615]
datasets = datasets.reindex(index=dataset_ids)


# In[65]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[66]:


data.head()


# ## Subset to the genes currently in SGD

# In[67]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[68]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[69]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[70]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[71]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[72]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[73]:


from IO.save_data_to_db3 import *


# In[74]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




