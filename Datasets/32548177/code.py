#!/usr/bin/env python
# coding: utf-8

# In[3]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[4]:


paper_pmid = 32548177
paper_name = 'edouarzin_vediyappan_2020' 


# In[5]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[6]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[8]:


original_data = pd.read_excel('raw_data/mic-07-146-s02.xls', sheet_name='12_27_11_15_06_18_nQuantile_nq', skiprows=2)


# In[9]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[10]:


original_data.head()


# In[11]:


original_data['orf'] = original_data['probeid'].astype(str)


# In[12]:


original_data['orf'] = original_data['orf'].apply(lambda x: x.split(':')[0])


# In[13]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[14]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[15]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[16]:


data_cols = ['(2) (Albicanol) EP6_0.05mg/ml','(1) (Drimenol) EP9_0.025mg/ml','(1) ( Drimenol) EP9_0.0375mg/ml']


# In[17]:


original_data.set_index('orf', inplace=True)


# In[18]:


original_data1 = original_data.loc[original_data['essential_gene']=='no', data_cols].copy()
original_data2 = original_data.loc[original_data['essential_gene']=='yes', data_cols].copy()


# In[19]:


original_data1 = original_data1.apply(pd.to_numeric, axis=1, errors='coerce')
original_data2 = original_data2.apply(pd.to_numeric, axis=1, errors='coerce')


# In[20]:


original_data1 = original_data1.groupby(original_data1.index).mean()
original_data2 = original_data2.groupby(original_data2.index).mean()


# In[21]:


original_data1.shape


# In[22]:


original_data2.shape


# In[24]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[25]:


original_data.head()


# # Prepare the final dataset

# In[26]:


data = original_data.copy()


# In[27]:


dataset_ids = [21872, 21865, 21870, 21871, 21866, 21869]
datasets = datasets.reindex(index=dataset_ids)


# In[28]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[29]:


data.head()


# ## Subset to the genes currently in SGD

# In[30]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[31]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[32]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[33]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[34]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[35]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[36]:


from IO.save_data_to_db3 import *


# In[37]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




