#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 16168084
paper_name = 'eide_harper_2005' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/13059_2005_1112_MOESM2_ESM.xlsx', sheet_name='Normalized concentrations')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data['orf'] = original_data['Gene Mutated'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[12]:


original_data.loc[original_data['orf']=='YILO16W','orf'] = 'YIL016W'


# In[13]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[14]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[15]:


original_data.set_index('orf', inplace=True)


# In[16]:


original_data.drop(columns='Gene Mutated', inplace=True)


# In[17]:


original_data = original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[18]:


original_data = original_data.groupby(original_data.index).mean()


# In[19]:


original_data.shape


# # Load dataset ids

# In[22]:


dt = pd.read_excel('extras/datasets.xlsx', sheet_name='Sheet1', header=None)
dt.head()


# In[23]:


dt.set_index(0, inplace=True)


# In[24]:


dt = dt.reindex(index=original_data.columns.values)


# In[26]:


dataset_ids = dt[1].values


# # Load & process tested strains

# In[27]:


tested = pd.read_excel('raw_data/13059_2005_1112_MOESM1_ESM.xlsx', sheet_name='Genes Analyzed')


# In[28]:


tested.head()


# In[29]:


tested['orf'] = tested['First pass (n=1) samples'].astype(str)


# In[30]:


tested['orf'] = clean_orf(tested['orf'])


# In[33]:


typo_fixes = {'YAL022C-2':'YAL022C','YAL022C-3':'YAL022C','YELOO1C':'YEL001C','YIL077C-3':'YIL077C','YOR134W-':'YOR134W'}


# In[34]:


tested['orf'] = tested['orf'].apply(lambda x: typo_fixes[x] if x in typo_fixes.keys() else x)


# In[35]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[36]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[37]:


tested = tested.loc[t,:]


# In[38]:


tested_orfs = tested['orf'].unique()


# In[39]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[40]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[41]:


data = original_data.copy()


# In[42]:


datasets = datasets.reindex(index=dataset_ids)


# In[43]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[44]:


data.head()


# ## Subset to the genes currently in SGD

# In[45]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[46]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

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

data_all.head()


# # Print out

# In[50]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[51]:


from IO.save_data_to_db3 import *


# In[52]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




