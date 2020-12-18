#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 32469861
paper_name = 'liu_liu_2020' 


# In[8]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[9]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[24]:


original_data = pd.read_excel('raw_data/journal.pgen.1008798.s007.xlsx', sheet_name='Combined data', skiprows=2)


# In[25]:


original_data.columns = ['rank1','orf1','','rank2','orf2','','rank3','orf3','']


# In[26]:


original_data.head()


# In[27]:


ranks = pd.concat([original_data['rank1'], original_data['rank2'], original_data['rank3']], axis=0, ignore_index=True)


# In[28]:


orfs = pd.concat([original_data['orf1'], original_data['orf2'], original_data['orf3']], axis=0, ignore_index=True)


# In[29]:


hit_data = ranks.to_frame().join(orfs.to_frame(), how='outer', lsuffix='_rank', rsuffix='_orf')


# In[30]:


hit_data.head()


# In[31]:


hit_data['0_orf'] = hit_data['0_orf'].astype(str)


# In[32]:


# Eliminate all white spaces & capitalize
hit_data['0_orf'] = clean_orf(hit_data['0_orf'])


# In[33]:


# Translate to ORFs 
hit_data['orfs'] = translate_sc(hit_data['0_orf'], to='orf')


# In[34]:


# Make sure everything translated ok
t = looks_like_orf(hit_data['orfs'])
print(hit_data.loc[~t,])


# In[35]:


hit_data = hit_data.loc[t,:]


# In[36]:


hit_data.set_index('orfs', inplace=True)
hit_data.index.name='orf'


# In[37]:


hit_data['data'] = 1


# In[38]:


hit_data = hit_data.groupby(hit_data.index).mean()


# In[39]:


hit_data.shape


# # Load & process tested strains

# In[40]:


tested_strains = pd.read_excel('raw_data/Original data after SGA Scoring sorted.xls', sheet_name='Combined data')


# In[41]:


tested_strains['Array ORF'] = tested_strains['Array ORF'].astype(str)


# In[42]:


# Eliminate all white spaces & capitalize
tested_strains['Array ORF'] = clean_orf(tested_strains['Array ORF'])


# In[43]:


# Translate to ORFs 
tested_strains['orfs'] = translate_sc(tested_strains['Array ORF'], to='orf')


# In[44]:


# Make sure everything translated ok
t = looks_like_orf(tested_strains['orfs'])
print(tested_strains.loc[~t,])


# In[45]:


tested = tested_strains['orfs'].unique()


# In[46]:


missing = [orf for orf in hit_data.index.values if orf not in tested]


# In[47]:


missing


# # Prepare the final dataset

# In[68]:


dataset_ids = [16543]


# In[69]:


datasets = datasets.reindex(index=dataset_ids)


# In[70]:


data = pd.DataFrame(index=tested, columns=datasets['name'].values, data=0)


# In[71]:


data.loc[hit_data.index, datasets['name'].values[0]] = hit_data['data']


# In[72]:


data = data.groupby(data.index).mean()


# In[73]:


# Create row index
data.index.name='orf'


# In[74]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[75]:


print('Final data dimensions: %d x %d' % (data.shape))


# ## Subset to the genes currently in SGD

# In[76]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[77]:


data.shape


# In[78]:


gene_ids.shape


# In[79]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[80]:


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


# In[84]:


data_all = data.join(data_norm)


# In[85]:


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


# In[89]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




