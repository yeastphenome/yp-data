#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 20936605
paper_name = 'yoshida_yoshimoto_2011' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[9]:


original_data1 = pd.read_csv('raw_data/hits_alpha.txt', header=None, names=['genes','data'], sep='\t')


# In[10]:


print('Original data dimensions: %d x %d' % (original_data1.shape))


# In[11]:


original_data1.head()


# In[13]:


original_data1['genes'] = original_data1['genes'].astype(str)


# In[15]:


# Eliminate all white spaces & capitalize
original_data1['genes'] = clean_genename(original_data1['genes'])


# In[20]:


original_data1.loc[original_data1['genes']=='ZSP1','genes'] = 'YBR287W'


# In[21]:


# Translate to ORFs 
original_data1['orf'] = translate_sc(original_data1['genes'], to='orf')


# In[22]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orf'])
print(original_data1.loc[~t,])


# In[23]:


original_data1['data'] = original_data1['data'].astype(int)


# In[24]:


original_data1.set_index('orf', inplace=True)


# In[25]:


original_data1 = original_data1[['data']].copy()


# In[26]:


original_data1 = original_data1.groupby(original_data1.index).mean()


# In[27]:


original_data1.shape


# # Load data (2)

# In[28]:


m = pd.read_csv('raw_data/hits_a.txt', header=None, sep='\t')


# In[30]:


m[0] = m[0].apply(lambda x: x.split(','))


# In[32]:


genes = [gene for r in m[0] for gene in r ]


# In[34]:


original_data2 = pd.DataFrame(data={'gene': genes})


# In[35]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[36]:


original_data2['gene'] = original_data2['gene'].astype(str)


# In[37]:


# Eliminate all white spaces & capitalize
original_data2['gene'] = clean_genename(original_data2['gene'])


# In[38]:


# Translate to ORFs 
original_data2['orf'] = translate_sc(original_data2['gene'], to='orf')


# In[39]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[40]:


original_data2['data'] = -1


# In[41]:


original_data2.set_index('orf', inplace=True)


# In[42]:


original_data2 = original_data2[['data']].copy()


# In[43]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[44]:


original_data2.shape


# # Load & process tested strains

# In[45]:


tested1 = pd.read_excel('raw_data/Mat_alpha_obs_v6.0.xls', sheet_name='DATA')


# In[46]:


tested1.head()


# In[47]:


tested1['orf'] = tested1['ORF name'].astype(str)


# In[48]:


tested1['orf'] = clean_orf(tested1['orf'])


# In[49]:


tested1['orf'] = translate_sc(tested1['orf'], to='orf')


# In[50]:


# Make sure everything translated ok
t = looks_like_orf(tested1['orf'])
print(tested1.loc[~t,])


# In[51]:


tested1 = tested1.loc[t,:]


# In[52]:


tested_orfs = tested1['orf'].unique()


# In[53]:


missing = [orf for orf in original_data1.index.values if orf not in tested_orfs]
missing


# In[54]:


tested_orfs = list(tested_orfs) + missing


# In[55]:


original_data1 = original_data1.reindex(index=tested_orfs, fill_value=0)


# # Load tested strains (2)

# In[ ]:


# Received 2 Mat-a files, will take the union of the strains


# In[68]:


files = ['Mat_a_obs_v5.0.xls','mat_a.xls']
sheets = ['DATA','Sheet1']


# In[69]:


tested_list = []
for ixf, f in enumerate(files):
    tested2 = pd.read_excel('raw_data/' + f, sheet_name=sheets[ixf])
    tested2['orf'] = tested2['ORF name'].astype(str)
    tested2['orf'] = clean_orf(tested2['orf'])
    tested2.loc[tested2['orf']=='YLR287-A','orf'] = 'YLR287C-A'
    tested2['orf'] = translate_sc(tested2['orf'], to='orf')
    # Make sure everything translated ok
    t = looks_like_orf(tested2['orf'])
    print(tested2.loc[~t,])
    tested2 = tested2.loc[t,:]
    tested_orfs2 = tested2['orf'].unique()
    tested_list.append(tested_orfs2)


# In[70]:


tested_orfs = list(tested_list[0]) + list(tested_list[1])


# In[71]:


tested_orfs = np.unique(np.array(tested_orfs))
tested_orfs.shape


# In[72]:


missing = [orf for orf in original_data2.index.values if orf not in tested_orfs]
missing


# In[73]:


tested_orfs = list(tested_orfs) + missing


# In[74]:


original_data2 = original_data2.reindex(index=tested_orfs, fill_value=0)


# # Merge

# In[75]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[76]:


original_data.head()


# # Prepare the final dataset

# In[77]:


data = original_data.copy()


# In[78]:


dataset_ids = [16694, 16695]
datasets = datasets.reindex(index=dataset_ids)


# In[79]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[80]:


data.head()


# ## Subset to the genes currently in SGD

# In[81]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[82]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[83]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[84]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[85]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[86]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[87]:


from IO.save_data_to_db3 import *


# In[88]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




