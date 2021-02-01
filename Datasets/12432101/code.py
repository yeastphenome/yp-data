#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 12432101
paper_name = 'deutschbauer_davis_2002' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[54]:


files = ['sporulation_deficient.txt','sporulation_proficient.txt','germination.txt']


# In[55]:


original_data_list = []
for f in files:
    original_data = pd.read_csv('raw_data/' + f, sep='\t')
    print('Original data dimensions: %d x %d' % (original_data.shape))
    print(original_data.head())
    original_data['orf'] = original_data['Orf'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data.set_index('orf', inplace=True)
    
    if f.startswith('sporulation'):
        original_data = original_data[['Spo1','Spo2']].apply(pd.to_numeric, axis=1, errors='coerce')
        original_data[['Spo1','Spo2']] = 1/original_data[['Spo1','Spo2']]
    else:
        original_data = original_data[['Germ1','Germ2']].apply(pd.to_numeric, axis=1, errors='coerce')
        
    original_data['data'] = original_data.mean(axis=1)
    
    original_data = original_data[['data']].copy()
    
    original_data = original_data.groupby(original_data.index).mean()
    
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[56]:


original_data1 = pd.concat(original_data_list[:2], axis=0)
original_data1 = original_data1.groupby(original_data1.index).mean()


# In[57]:


original_data2 = original_data_list[2]


# In[58]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[59]:


original_data.shape


# In[60]:


original_data.head()


# # Load & process tested strains

# In[61]:


tested = pd.read_csv('raw_data/SpoGerm_RawData.txt', sep='\t')


# In[62]:


tested.head()


# In[63]:


tested['orf'] = tested['ORF     '].astype(str)


# In[64]:


tested['orf'] = clean_orf(tested['orf'])


# In[65]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[66]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[67]:


tested['essential'] = is_essential(tested['orf'])


# In[68]:


tested = tested.loc[~tested['essential']]


# In[69]:


tested.shape


# In[70]:


tested_orfs = tested['orf'].unique()


# In[71]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[72]:


# WT = 1
original_data = original_data.reindex(index=tested_orfs, fill_value=1)


# In[73]:


original_data.shape


# # Prepare the final dataset

# In[74]:


data = original_data.copy()


# In[75]:


dataset_ids = [478, 16004]
datasets = datasets.reindex(index=dataset_ids)


# In[76]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[77]:


data.head()


# ## Subset to the genes currently in SGD

# In[78]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[79]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[80]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[81]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[82]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[83]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[52]:


from IO.save_data_to_db3 import *


# In[53]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




