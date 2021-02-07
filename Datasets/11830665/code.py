#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 11830665
paper_name = 'fleming_blackman_2002' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


files = ['hit_data.txt','hit_data2.txt']


# In[64]:


original_data_list = []
for ixf, f in enumerate(files):
    original_data = pd.read_csv('raw_data/' + f, sep='\t')
    print('Original data dimensions: %d x %d' % (original_data.shape))
    print(original_data.head())
    original_data['gene'] = original_data['genenames'].astype(str)
    original_data['gene'] = clean_genename(original_data['gene'])
    original_data['orf'] = translate_sc(original_data['gene'], to='orf')
    # Make sure everything translated ok
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    
    original_data.set_index('orf', inplace=True)
    
    if ixf == 1:
        data_cols = ['ps-519','ps-341']
        for c in data_cols:
            original_data[c] = 1
    else:
        data_cols = ['ps-519','ps-341']
        original_data[data_cols] = original_data[data_cols].apply(pd.to_numeric, axis=1, errors='coerce')
        
        wt = original_data.loc['WT',data_cols]
        
        original_data[data_cols] = original_data[data_cols].subtract(wt, axis=1).divide(wt, axis=1)
        original_data.drop(index='WT', inplace=True)
        
    original_data = original_data[data_cols].copy()
    original_data = original_data.groupby(original_data.index).mean()
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[65]:


original_data = pd.concat(original_data_list, axis=0)


# In[66]:


original_data.shape


# In[67]:


original_data = original_data.groupby(original_data.index).mean()


# In[68]:


original_data.shape


# # Load & process tested strains

# In[69]:


tested = pd.read_excel('raw_data/Table 2 (web supplement).xlsx', sheet_name='Sheet1', skiprows=4)


# In[70]:


tested.head()


# In[71]:


tested.shape


# In[72]:


tested = tested.loc[tested['Genomics']=='+',:]
tested.shape


# In[73]:


tested['orf'] = tested['ORF Name'].astype(str)


# In[74]:


tested['orf'] = clean_orf(tested['orf'])


# In[75]:


# Insert missing '-'
tested['orf'] = tested['orf'].apply(lambda x: x[:7] + '-' + x[7] if len(x) == 8 else x)


# In[76]:


tested['orf'] = translate_sc(tested['orf'].values, to='orf')


# In[77]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[78]:


tested_orfs = tested['orf'].unique()


# In[79]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[80]:


tested_orfs = list(tested_orfs) + missing


# In[81]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[82]:


data = original_data.copy()


# In[83]:


dataset_ids = [1308, 471]
datasets = datasets.reindex(index=dataset_ids)


# In[84]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[85]:


data.head()


# ## Subset to the genes currently in SGD

# In[86]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[87]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[88]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[89]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[90]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[91]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[92]:


from IO.save_data_to_db3 import *


# In[93]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




