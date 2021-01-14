#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 19853568
paper_name = 'carroll_drubin_2009' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/1-s2.0-S1534580709003438-mmc2-2.xls', sheen_name='Table S1', skiprows=7, header=None)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data = original_data.loc[original_data[1]=='KO']


# In[8]:


original_data.head()


# In[9]:


original_data['ORF'] = original_data[3].astype(str)


# In[10]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[11]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'].values, to='orf')


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[13]:


original_data[4] = [s.strip() for s in original_data[4]]


# In[15]:


data_vals_dict = {'++': 2, '+': 1, '-': -1, '- -': -2}


# In[16]:


original_data['data'] = original_data[4].apply(lambda x: data_vals_dict[x])


# In[17]:


original_data.set_index('ORF', inplace=True)
original_data.index.name='orf'


# In[18]:


original_data = original_data[['data']].copy()


# In[19]:


original_data = original_data.groupby(original_data.index).mean()


# In[20]:


original_data.shape


# # Prepare the final dataset

# In[21]:


data = original_data.copy()


# In[22]:


dataset_ids = [16455]
datasets = datasets.reindex(index=dataset_ids)


# In[23]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[24]:


data.head()


# ## Subset to the genes currently in SGD

# In[25]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[26]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[27]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[28]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[29]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[30]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[31]:


from IO.save_data_to_db3 import *


# In[32]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




