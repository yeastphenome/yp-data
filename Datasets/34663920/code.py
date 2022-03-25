#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 34663920
paper_name = 'vieitez_beltrao_2021' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/41587_2021_1051_MOESM5_ESM.xlsx', sheet_name='Table S3 – S_Scores of chemical')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data = original_data.loc[original_data['Mutant_type']=='KO']


# In[9]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[10]:


t = pd.pivot_table(original_data, index='Systematic_name', columns='Condition', values='Score')


# In[14]:


original_data = t.copy().reset_index()


# In[15]:


original_data['orf'] = original_data['Systematic_name'].astype(str)


# In[16]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[17]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[18]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[19]:


original_data.loc[original_data['orf']=='YLR287-A','orf'] = 'YLR287C-A'


# In[20]:


original_data.set_index('orf', inplace=True)


# In[21]:


original_data = original_data.groupby(original_data.index).mean()


# In[22]:


original_data.shape


# In[23]:


original_data.head()


# # Prepare the final dataset

# In[52]:


data = original_data.copy()


# In[58]:


cond2dt = pd.read_csv('raw_data/condition_2_dataset.txt', sep='\t')
print(cond2dt.shape)


# In[59]:


cond2dt['condition_id_fix'] = cond2dt['condition_id_fix'].str.upper()


# In[60]:


cond2dt.set_index('condition_id_fix', inplace=True)


# In[61]:


cond2dt = cond2dt.reindex(index=data.columns)


# In[63]:


cond2dt.loc[cond2dt['dataset_id'].isnull()]


# In[64]:


cond2dt.head()


# In[65]:


dataset_ids = cond2dt['dataset_id'].values
datasets = datasets.reindex(index=dataset_ids)


# In[66]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[67]:


data.head()


# ## Subset to the genes currently in SGD

# In[68]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[69]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[70]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[71]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[72]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[73]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[74]:


from IO.save_data_to_db3 import *


# In[75]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




