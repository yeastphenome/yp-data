#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 16507139
paper_name = 'narayanaswamy_marcotte_2006' 


# In[86]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[87]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[89]:


original_data = pd.read_excel('raw_data/gb-2006-7-1-r6-s2.xlsx',sheet_name='Cellma Grades', skiprows=39)


# In[90]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[91]:


original_data.head()


# In[92]:


original_data['orf'] = original_data['GENE'].astype(str)


# In[93]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[94]:


original_data.loc[original_data['orf']=='YLR287-A','orf'] = 'YLR287C-A'


# In[95]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[96]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[97]:


original_data.set_index('orf', inplace=True)
original_data.drop(columns=['GENE'], inplace=True)


# In[98]:


original_data = original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[99]:


original_data = original_data.groupby(original_data.index).mean()


# In[100]:


original_data.shape


# In[101]:


# Multiply intensity by penetrance
int_cols = [c for c in original_data.columns.values if '_INT' in c]
pen_cols = [c for c in original_data.columns.values if '_PEN' in c]


# In[102]:


for ic,c in enumerate(int_cols):
    original_data[int_cols[ic]] = original_data[int_cols[ic]] * original_data[pen_cols[ic]] / 4


# In[103]:


original_data.drop(columns=pen_cols, inplace=True)


# In[104]:


# Average between observers
obs1_cols = [c for c in original_data.columns.values if '22_' in c]
obs2_cols = [c for c in original_data.columns.values if '24_' in c]


# In[105]:


for ic, c in enumerate(obs1_cols):
    original_data[obs1_cols[ic]] = (original_data[obs1_cols[ic]] + original_data[obs2_cols[ic]]) / 2


# In[106]:


original_data.drop(columns=obs2_cols, inplace=True)


# # Load dataset_ids

# In[107]:


mp = pd.read_excel('raw_data/phenotype_mapping.xlsx', sheet_name='Sheet1', skiprows=1)


# In[108]:


mp.set_index('Unnamed: 0', inplace=True)


# In[111]:


mp.columns = [int(c) for c in mp.columns]


# In[112]:


original_data = original_data.reindex(columns=mp.index.values)


# In[113]:


original_data.head()


# In[114]:


original_data2 = pd.DataFrame(index=original_data.index, columns=mp.columns, data=0)


# In[115]:


for p in mp.index.values:
    for d in mp.columns.values:
        v = mp.loc[p,d]
        if np.abs(v) > 0:
            original_data2.loc[:,d] = original_data2.loc[:,d] + original_data.loc[:,p] * v


# In[116]:


original_data2.head()


# # Prepare the final dataset

# In[117]:


data = original_data2.copy()


# In[118]:


dataset_ids = data.columns.values
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


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[127]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[128]:


from IO.save_data_to_db3 import *


# In[129]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




