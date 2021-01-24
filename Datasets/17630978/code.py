#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 17630978
paper_name = 'pagani_arino_2007' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, names=['genes','data'], sep='\t')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data['genes'] = original_data['genes'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[9]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[11]:


original_data.set_index('orfs', inplace=True)


# In[12]:


original_data.index.name = 'orf'


# In[14]:


original_data['data'] = pd.to_numeric(original_data['data'], errors='coerce')


# In[15]:


original_data = original_data[['data']].copy()


# In[16]:


original_data = original_data.groupby(original_data.index).mean()


# In[17]:


original_data.shape


# # Load & process tested strains

# In[18]:


tested = pd.read_excel('raw_data/EUROFAN haploid collection .xlsx', sheet_name='GENERAL 07-02-11')


# In[19]:


tested['orf'] = tested['Systematic Name '].astype(str)
tested['orf'] = clean_orf(tested['orf'])


# In[20]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[21]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[22]:


tested.loc[~t,'orf'].unique()


# In[23]:


tested = tested.loc[t,:]


# In[24]:


tested_orfs = tested['orf'].unique()


# In[25]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[26]:


# Decided to add the 1 missing strain
tested_orfs = list(tested_orfs) + missing


# In[27]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[28]:


data = original_data.copy()


# In[29]:


dataset_ids = [16619]
datasets = datasets.reindex(index=dataset_ids)


# In[30]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[31]:


data.head()


# ## Subset to the genes currently in SGD

# In[32]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[33]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[34]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[35]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[36]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[37]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[38]:


from IO.save_data_to_db3 import *


# In[39]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




