#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 21960436
paper_name = 'dos_santos_sa_correia_2011' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[9]:


files = ['hits_genenames_hs.txt','hits_genenames_s.txt','hits_genenames_r.txt']
scores = [-2,-1,1]


# In[10]:


original_data_list = []
for ixf, f in enumerate(files):
    original_data = pd.read_csv('raw_data/' + f, header=None, names=['gene'])
    print('Original data dimensions: %d x %d' % (original_data.shape))
    original_data['gene'] = original_data['gene'].astype(str)
    # Eliminate all white spaces & capitalize
    original_data['gene'] = clean_orf(original_data['gene'])
    # Translate to ORFs 
    original_data['orf'] = translate_sc(original_data['gene'], to='orf')
    
    # Make sure everything translated ok
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    
    original_data['data'] = scores[ixf]
    original_data.set_index('orf', inplace=True)
    original_data = original_data[['data']].copy()
    original_data = original_data.groupby(original_data.index).mean()
    
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[11]:


original_data = pd.concat(original_data_list, axis=1)


# In[13]:


original_data.columns = ['HS','S','R']


# In[14]:


original_data.head()


# In[9]:


# Adjust the overlapping strains as follows: eliminate HS from S; eliminate
# HS ^ R and S ^ R from both


# In[17]:


original_data.loc[original_data[['HS','R']].notnull().sum(axis=1) > 1,['HS','R']] = np.nan


# In[18]:


original_data.loc[original_data[['S','R']].notnull().sum(axis=1) > 1,['S','R']] = np.nan


# In[19]:


original_data.loc[original_data[['HS','S']].notnull().sum(axis=1) > 1,['S']] = np.nan


# In[21]:


original_data['data'] = original_data.sum(axis=1)


# In[23]:


original_data = original_data[['data']].copy()


# In[24]:


original_data = original_data.groupby(original_data.index).mean()


# In[25]:


original_data.shape


# # Load & process tested strains

# In[26]:


tested = pd.read_excel('raw_data/List of strains tested.xlsx', sheet_name='Tabelle2')


# In[27]:


tested.head()


# In[28]:


tested['orf'] = tested['ORF'].astype(str)


# In[29]:


tested['orf'] = clean_orf(tested['orf'])


# In[30]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[31]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[32]:


tested_orfs = tested['orf'].unique()


# In[33]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[34]:


tested_orfs = list(tested_orfs) + missing


# In[35]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# In[42]:


original_data.index.name='orf'


# # Prepare the final dataset

# In[43]:


data = original_data.copy()


# In[44]:


dataset_ids = [149]
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




