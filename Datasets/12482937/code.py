#!/usr/bin/env python
# coding: utf-8

# In[14]:


get_ipython().run_line_magic('run', '../yp_utils.py')
import re


# # Initial setup

# In[2]:


paper_pmid = 12482937
paper_name = 'chang_brown_2002' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[7]:


original_data = pd.read_csv('raw_data/MMS Sensitive Strains.txt', header=None, sep='\n')


# In[8]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[9]:


original_data.head()


# In[18]:


original_data['gene'] = original_data[0].apply(lambda x: re.split(' |\/',x)[0])
original_data['data'] = original_data[0].apply(lambda x: re.split(' ',x)[-1])


# In[19]:


original_data.head()


# In[20]:


original_data['gene'] = original_data['gene'].astype(str)


# In[21]:


# Eliminate all white spaces & capitalize
original_data['gene'] = clean_genename(original_data['gene'])


# In[24]:


typo_fixes = {'KIM3': 'MMS1', 'SW16': 'SWI6'}
original_data['gene'] = original_data['gene'].apply(lambda x: typo_fixes[x] if x in typo_fixes.keys() else x)


# In[25]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['gene'], to='orf')


# In[26]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[27]:


original_data['data'] = original_data['data'].apply(lambda x: -len(x.strip()))


# In[28]:


original_data.set_index('orf', inplace=True)


# In[29]:


original_data = original_data[['data']].copy()


# In[30]:


original_data = original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[31]:


original_data = original_data.groupby(original_data.index).mean()


# In[32]:


original_data.shape


# # Load & process tested strains

# In[33]:


tested = pd.read_excel('raw_data/Old array position-before Dec 8, 05.xlsx', sheet_name='Sheet1')


# In[34]:


tested.head()


# In[35]:


tested['orf'] = tested['Systematic Name'].astype(str)


# In[36]:


tested['orf'] = clean_orf(tested['orf'])


# In[39]:


tested.loc[tested['orf']=='YPL072WA','orf'] = 'YPL072W-A'


# In[40]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[41]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[42]:


tested = tested.loc[t,:]


# In[43]:


tested_orfs = tested['orf'].unique()


# In[44]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[45]:


tested_orfs = list(tested_orfs) + missing


# In[46]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[47]:


data = original_data.copy()


# In[48]:


dataset_ids = [65]
datasets = datasets.reindex(index=dataset_ids)


# In[49]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[50]:


data.head()


# ## Subset to the genes currently in SGD

# In[51]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[52]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[53]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[54]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[55]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[56]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[57]:


from IO.save_data_to_db3 import *


# In[58]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




