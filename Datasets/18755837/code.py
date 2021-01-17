#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 18755837
paper_name = 'huang_bystrom_2008' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[53]:


original_data = pd.read_csv('raw_data/Table1-2.txt', header=None, names=['genes','data'], sep='\t')


# In[54]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[55]:


original_data['genes'] = original_data['genes'].astype(str)


# In[56]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[57]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[58]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[59]:


original_data = original_data.loc[t,:]


# In[60]:


original_data.set_index('orfs', inplace=True)


# In[61]:


original_data.index.name='orf'


# In[62]:


original_data['data'] = original_data['data'].astype(int)


# In[63]:


original_data = original_data[['data']].copy()


# # Load & process tested strains

# In[64]:


orf_pattern = 'Y[A-P][RL][0-9]{3}[CW](-[A-H])*'
p = re.compile(orf_pattern)


# In[65]:


file_path = 'raw_data/Deletion collection homo dipl.txt'

tested_orfs = []
with open(file_path,'r') as fp:
    line = fp.readline()
    while line:
        res = p.search(line)
        if res:
            tested_orfs.append(res.group(0))
        line = fp.readline()
        


# In[66]:


tested = pd.DataFrame(data={'orf': tested_orfs})


# In[67]:


tested['orf'] = clean_orf(tested['orf'])


# In[68]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[69]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[70]:


tested_orfs = tested['orf'].unique()


# In[71]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]


# In[72]:


missing


# In[73]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# In[74]:


original_data.head()


# # Prepare the final dataset

# In[75]:


data = original_data.copy()


# In[76]:


dataset_ids = [16436]
datasets = datasets.reindex(index=dataset_ids)


# In[77]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[78]:


data.head()


# ## Subset to the genes currently in SGD

# In[79]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[80]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[81]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[82]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[83]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[85]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[86]:


from IO.save_data_to_db3 import *


# In[87]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




