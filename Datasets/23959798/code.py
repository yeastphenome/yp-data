#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 23959798
paper_name = 'dong_rutherford_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/zmb999100152sd1.xlsx', sheet_name='Sheet1', skiprows=1)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


all_genes = []
for gene_list in original_data['Gene(s)'].values:
    all_genes = all_genes + gene_list.split(',')


# In[9]:


all_genes = list(set(all_genes))


# In[10]:


len(all_genes)


# In[11]:


# Eliminate all white spaces & capitalize
all_genes = clean_genename(all_genes)


# In[12]:


# Translate to ORFs 
all_orfs = translate_sc(all_genes, to='orf')


# In[13]:


# Make sure everything translated ok
t = looks_like_orf(all_orfs)
print(np.array(all_orfs)[~np.array(t)])


# # Load & process tested strains

# In[14]:


tested = pd.read_excel('raw_data/Gene list.xls', sheet_name='Sheet1', header=None)


# In[15]:


tested[1] = clean_orf(tested[1])


# In[16]:


tested[1] = translate_sc(tested[1], to='orf')


# In[17]:


# Make sure everything translated ok
t = looks_like_orf(tested[1])
print(tested.loc[~t])


# In[18]:


tested_orfs = tested[1].unique()


# In[22]:


missing = [orf for orf in all_orfs if orf not in tested_orfs]


# In[23]:


missing


# # Prepare the final dataset

# In[24]:


dataset_ids = [16559]


# In[25]:


datasets = datasets.reindex(index=dataset_ids)


# In[26]:


data = pd.DataFrame(index=tested_orfs, columns=dataset_ids, data=0)


# In[27]:


data.loc[all_orfs] = 1


# In[28]:


data = data.groupby(data.index).mean()


# In[29]:


# Create row index
data.index.name='orf'


# In[32]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[33]:


data.sum(axis=0)


# ## Subset to the genes currently in SGD

# In[34]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[35]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[36]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[37]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[38]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[39]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[40]:


from IO.save_data_to_db3 import *


# In[41]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




