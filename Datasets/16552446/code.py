#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 16552446
paper_name = 'gatbonton_bedalov_2006' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[84]:


original_data = pd.read_excel('raw_data/gatbonton_bedalov_2006_hits.xlsx', sheet_name='Sheet1', header=None)


# In[85]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[86]:


original_data.head()


# In[87]:


original_data['gene'] = original_data[0].astype(str)


# In[88]:


# Eliminate all white spaces & capitalize
original_data['gene'] = clean_genename(original_data['gene'])


# In[89]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['gene'], to='orf')


# In[90]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[91]:


original_data = original_data.loc[t,]


# In[124]:


# Converting the scores such that:
# a) 1 = weakest phenotype, 3 = strongest phenotype
# b) long telomeres = positive scores, short telomeres = negative scores


# In[92]:


original_data['data'] = original_data[2].astype(int)-4


# In[96]:


original_data.loc[original_data[1] == 'L','data'] = -original_data.loc[original_data[1] == 'L','data']


# In[98]:


original_data.set_index('orf', inplace=True)


# In[99]:


original_data = original_data[['data']].copy()


# In[100]:


original_data = original_data.groupby(original_data.index).mean()


# In[101]:


original_data.shape


# # Load & process tested strains

# In[102]:


tested = pd.read_excel('raw_data/genelist_altered_020806.xlsx', sheet_name='mat alpha copy.txt')


# In[103]:


tested.head()


# In[104]:


tested['orf'] = tested['ORF name'].astype(str)


# In[105]:


tested['orf'] = clean_orf(tested['orf'])


# In[106]:


tested.loc[tested['orf']=='YYKL138C','orf'] = 'YKL138C'


# In[107]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[108]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[109]:


tested_orfs = tested['orf'].unique()


# In[110]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[111]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[112]:


data = original_data.copy()


# In[113]:


dataset_ids = [104]
datasets = datasets.reindex(index=dataset_ids)


# In[114]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[115]:


data.head()


# ## Subset to the genes currently in SGD

# In[116]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[117]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[118]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[119]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[120]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[121]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[122]:


from IO.save_data_to_db3 import *


# In[123]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




