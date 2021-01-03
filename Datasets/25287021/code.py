#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 25287021
paper_name = 'pereira_domingues_2014' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[7]:


files = ['10295_2014_1519_MOESM1_ESM.xlsx','10295_2014_1519_MOESM2_ESM.xlsx']
sheets = ['Table S1','Table S2']


# In[12]:


original_data_list = []
for ixf, f in enumerate(files):
    original_data = pd.read_excel('raw_data/' + f, sheet_name=sheets[ixf], skiprows=9, header=None)
    print('Original data dimensions: %d x %d' % (original_data.shape))
    original_data['genes'] = original_data[0].astype(str)
    original_data['genes'] = clean_genename(original_data['genes'])
    original_data['orf'] = translate_sc(original_data['genes'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    
    original_data = original_data.loc[t,:]
    
    original_data.loc[original_data[10].isnull(),10] = ''
    original_data['data'] = original_data[10].apply(lambda x: -len(x))
    
    original_data.set_index('orf', inplace=True)
    original_data = original_data[['data']].copy()
    
    original_data = original_data.groupby(original_data.index).mean()
    
    original_data_list.append(original_data)


# In[59]:


original_data = pd.concat(original_data_list, axis=1)


# In[60]:


original_data.index.name = 'orf'


# In[61]:


original_data[original_data.isnull()] = 0


# In[62]:


original_data.head()


# In[63]:


original_data.shape


# # Load & process tested strains

# In[64]:


tested = pd.read_excel('raw_data/chemogenomics.xlsx', sheet_name='hom.z_tdist_pval_nm.smallmol.co')


# In[65]:


tested['orf'] = tested['Orf '].astype(str)


# In[66]:


tested['orf'] = clean_orf(tested['orf'])


# In[67]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[68]:


t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[69]:


tested_orfs = np.unique(tested['orf'].values)


# In[70]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[71]:


# Add the missing strains
tested_orfs = list(tested_orfs) + missing


# In[72]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[73]:


data = original_data.copy()


# In[74]:


dataset_ids = [751, 752]
datasets = datasets.reindex(index=dataset_ids)


# In[75]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[76]:


data.head()


# ## Subset to the genes currently in SGD

# In[77]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[78]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[79]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[80]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[81]:


data_norm[data.isnull()] = np.nan
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


save_data_to_db(data_all, paper_pmid)


# In[ ]:




