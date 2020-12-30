#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 26466677
paper_name = 'stauffer_powers_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[9]:


original_data = pd.read_csv('raw_data/Table1.txt', sep='\n', skiprows=1)


# In[10]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[11]:


original_data.head()


# In[12]:


original_data['genes'] = original_data['Vacuolar\tProteins'].astype(str)


# In[13]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[16]:


original_data.loc[original_data['genes']=='RPL19-A','genes'] = 'RPL19A'


# In[17]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['genes'], to='orf')


# In[18]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[19]:


original_data = original_data.loc[t,:]


# In[20]:


original_data['data'] = -1


# In[21]:


original_data.set_index('orf', inplace=True)


# In[22]:


original_data = original_data[['data']].copy()


# In[23]:


original_data = original_data.groupby(original_data.index).mean()


# In[24]:


original_data.shape


# # Load & process tested strains

# In[25]:


tested = pd.read_excel('raw_data/Deletion Collection List for Anastasia.xlsx', sheet_name='Sheet1', skiprows=2)


# In[26]:


tested.head()


# In[33]:


tested['orf'] = tested['ORF'].astype(str)


# In[34]:


tested['orf'] = clean_orf(tested['orf'])


# In[35]:


tested.loc[tested['orf']=='YLR287-A','orf'] = 'YLR287C-A'


# In[36]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[37]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[38]:


tested.loc[~t,'orf'].values


# In[39]:


tested = tested.loc[t,:]


# In[40]:


tested_orfs = np.unique(tested['orf'].values)


# In[41]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[42]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[43]:


data = original_data.copy()


# In[44]:


dataset_ids = [766]
datasets = datasets.reindex(index=dataset_ids)


# In[45]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[46]:


data.head()


# ## Subset to the genes currently in SGD

# In[47]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[48]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[49]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[50]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[51]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[52]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[53]:


from IO.save_data_to_db3 import *


# In[54]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




