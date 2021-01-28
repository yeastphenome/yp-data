#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 15908144
paper_name = 'luban_schmidt_2005' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[73]:


original_data = pd.read_csv('raw_data/list_of_pet_mutants.txt', sep='\t', header=None)


# In[74]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[75]:


original_data.head()


# In[76]:


original_data['orf'] = original_data[0].astype(str)


# In[77]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[78]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[79]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[80]:


original_data['data'] = -1


# In[81]:


original_data.set_index('orf', inplace=True)


# In[82]:


original_data = original_data[['data']].copy()


# In[83]:


original_data = original_data.groupby(original_data.index).mean()


# In[84]:


original_data.shape


# # Load & process tested strains

# In[64]:


tested = pd.read_csv('raw_data/list_of_used_knockouts_PhD_Thesis_Luban.txt', sep='\t', header=None)


# In[65]:


tested.head()


# In[66]:


tested['orf'] = tested[0].astype(str)


# In[67]:


# Special cleanup
tested['orf'] = tested['orf'].apply(lambda x: x.replace('Yl','YI'))


# In[68]:


tested['orf'] = clean_orf(tested['orf'])


# In[69]:


typo_fixes = {'YBL098V':'YBL098W','YDR07SW':'YDR075W','YDR27SW':'YDR275W',
              'YDR51SW':'YDR515W','YDRS41C':'YDR541C','YEL0I6C':'YEL016C',
              'YIIL016C':'YHL016C','YIIL017W':'YHL017W','YHL0I9C':'YHL019C',
              'YJR09JC':'YJR091C','YNL09SC':'YNL095C','YPLOI8W':'YPL018W','YPL07LC':'YPL071C'}

tested['orf'] = tested['orf'].apply(lambda x: typo_fixes[x] if x in typo_fixes.keys() else x)


# In[70]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[71]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[72]:


tested_orfs = tested['orf'].unique()


# In[85]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[86]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[87]:


data = original_data.copy()


# In[88]:


dataset_ids = [417]
datasets = datasets.reindex(index=dataset_ids)


# In[89]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[90]:


data.head()


# ## Subset to the genes currently in SGD

# In[91]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[94]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[95]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[96]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[97]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[98]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[99]:


from IO.save_data_to_db3 import *


# In[100]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




