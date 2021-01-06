#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 23445507
paper_name = 'troppens_morrissey_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/Table S1.xlsx', sheet_name='Sheet1')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data['orf'] = original_data['Systematic name'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[9]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[10]:


original_data.loc[original_data['orf']=='YDR369', 'orf'] = 'YDR369C'
original_data.loc[original_data['orf']=='YOR382', 'orf'] = 'YOR382W'


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[12]:


data_dict = {'+': 1,'++': 2,'-': -1, 'nd': np.nan,'+/-': 0, '(+)/-': 0}


# In[13]:


original_data['data'] = original_data['120 μg/ml DAPG']


# In[14]:


original_data['data'] = [data_dict[x] for x in original_data['data']]


# In[15]:


original_data.set_index('orf', inplace=True)


# In[16]:


original_data = original_data[['data']].copy()


# In[17]:


original_data = original_data.groupby(original_data.index).mean()


# In[18]:


original_data.shape


# # Prepare the final dataset

# In[19]:


data = original_data.copy()


# In[20]:


dataset_ids = [16562]
datasets = datasets.reindex(index=dataset_ids)


# In[21]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[22]:


data.head()


# ## Subset to the genes currently in SGD

# In[23]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[24]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[25]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[26]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[27]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[28]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[29]:


from IO.save_data_to_db3 import *


# In[30]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




