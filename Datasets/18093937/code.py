#!/usr/bin/env python
# coding: utf-8

# In[82]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[83]:


paper_pmid = 18093937
paper_name = 'szymanski_goodman_2007' 


# In[84]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[85]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[86]:


original_data = pd.read_excel('raw_data//Table2.xlsx', sheet_name='Sheet1')


# In[87]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[88]:


original_data.head()


# In[89]:


original_data['gene'] = original_data['Unnamed: 0'].astype(str)


# In[90]:


# Eliminate all white spaces & capitalize
original_data['gene'] = clean_genename(original_data['gene'])


# In[91]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['gene'], to='orf')


# In[92]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[93]:


ph = pd.read_excel('raw_data/phenotype_mapping.xlsx', sheet_name='Sheet1')
ph.set_index('Phenotypes', inplace=True)
ph.head()


# In[94]:


original_data.set_index('orf', inplace=True)


# In[95]:


original_data.head()


# In[96]:


def ph_to_vec(s):
    p1 = [x.lower().strip() for x in s.split(',')]
    ph1 = ph.loc[p1,:].sum(axis=0)
    return ph1.to_list()


# In[97]:


original_data['data1'] = original_data['Log phase'].apply(lambda x: ph_to_vec(x))


# In[98]:


original_data['data2'] = original_data['Stationary phase'].apply(lambda x: ph_to_vec(x))


# In[99]:


original_data1 = pd.DataFrame(original_data['data1'].to_list(), index=original_data.index, columns=ph.columns)
original_data2 = pd.DataFrame(original_data['data2'].to_list(), index=original_data.index, columns=ph.columns)


# In[100]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[101]:


original_data = original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[102]:


original_data = original_data.groupby(original_data.index).mean()


# In[103]:


original_data.shape


# In[104]:


# Remove the phenotypes/datasets that (as a result of the conversion from categorical to binary phenotypes) ended up not having any hits
ix_cols = np.ravel(np.argwhere((np.abs(original_data) > 0).sum(axis=0).values == 0))
original_data.drop(columns=original_data.columns[ix_cols], inplace=True)


# In[105]:


original_data.shape


# # Prepare the final dataset

# In[106]:


data = original_data.copy()


# In[107]:


dataset_ids = [170, 701, 703, 705, 706, 707, 5387, 5388, 5389, 5390, 5391, 5392, 5393]
datasets = datasets.reindex(index=dataset_ids)


# In[108]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[109]:


data.head()


# ## Subset to the genes currently in SGD

# In[110]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[111]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[112]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[113]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[114]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[115]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[116]:


from IO.save_data_to_db3 import *


# In[117]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




