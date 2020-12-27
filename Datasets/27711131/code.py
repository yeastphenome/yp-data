#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 27711131
paper_name = 'jakubkova_tomaska_2016' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[7]:


original_data = pd.read_excel('raw_data/journal.pone.0164175.s006.xlsx', sheet_name='Tab S2', skiprows=1)


# In[8]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[9]:


original_data.head()


# In[10]:


original_data['orf'] = original_data['ORF'].astype(str)


# In[11]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[12]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[13]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[14]:


original_data = original_data.loc[t,:]


# In[16]:


original_data.set_index('orf', inplace=True)


# In[20]:


original_data['Nig'] = 0
original_data['Val'] = 0
p = 'phenotype of the deletant on Val/Nig'


# In[21]:


original_data.loc[original_data[p]=='Nig hypersensitive','Nig'] = -1
original_data.loc[original_data[p]=='Nig resistant','Nig'] = 1
original_data.loc[original_data[p]=='Val and Nig hypersensitive',:] = -1
original_data.loc[original_data[p]=='Val and Nig resistant',:] = 1
original_data.loc[original_data[p]=='Val hypersensitive','Val'] = -1
original_data.loc[original_data[p]=='Val resistant','Val'] = 1


# In[24]:


original_data = original_data[['Nig','Val']].copy()


# In[25]:


original_data = original_data.groupby(original_data.index).mean()


# In[26]:


original_data.shape


# # Load & process tested strains

# In[29]:


tested = pd.read_csv('raw_data/strain_a_mating_type.txt', sep='\t', skiprows=2, header=None)


# In[30]:


tested.head()


# In[31]:


tested['orf'] = tested[1].astype(str)


# In[32]:


tested['orf'] = clean_orf(tested['orf'])


# In[33]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[34]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[35]:


tested_orfs = np.unique(tested['orf'].values)


# In[37]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[38]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# In[39]:


original_data.shape


# # Prepare the final dataset

# In[40]:


data = original_data.copy()


# In[41]:


dataset_ids = [5182, 5176]
datasets = datasets.reindex(index=dataset_ids)


# In[42]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[43]:


data.head()


# ## Subset to the genes currently in SGD

# In[44]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[45]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[46]:


data.head()


# # Normalize

# In[47]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[48]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[49]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[50]:


data_all.head()


# # Print out

# In[51]:


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




