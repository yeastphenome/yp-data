#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 27146641
paper_name = 'johnson_wu_2016' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/c6mt00039h1.xlsx', sheet_name='Sheet1', skiprows=1, 
                              names=['orf','h2o_t0','h2o_t16','chr5_t0','chr5_t16','chr1_t0','chr1_t16'])


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data['orf'] = original_data['orf'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[10]:


# If possible, fix typos, omissions, etc.
original_data.loc[original_data['orf'].str.contains('BY4743AVERAGEN128'),'orf'] = 'WT'


# In[11]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])


# In[13]:


print(original_data.loc[~t,])


# In[14]:


# Normalize by t16 by t0, treated vs untreated and mut vs wt
original_data['h2o_ratio'] = original_data['h2o_t16'] / original_data['h2o_t0']
original_data['chr5_ratio'] = original_data['chr5_t16'] / original_data['chr5_t0']
original_data['chr1_ratio'] = original_data['chr1_t16'] / original_data['chr1_t0']


# In[15]:


original_data['h2o_ratio_wt'] = original_data['h2o_ratio'] / original_data.loc[original_data['orf']=='WT','h2o_ratio'].values
original_data['chr5_ratio_wt'] = original_data['chr5_ratio'] / original_data.loc[original_data['orf']=='WT','chr5_ratio'].values
original_data['chr1_ratio_wt'] = original_data['chr1_ratio'] / original_data.loc[original_data['orf']=='WT','chr1_ratio'].values


# In[16]:


original_data['chr5_ratio_wt_unt'] = original_data['chr5_ratio_wt'] / original_data['h2o_ratio_wt']
original_data['chr1_ratio_wt_unt'] = original_data['chr1_ratio_wt'] / original_data['h2o_ratio_wt']


# In[17]:


original_data.head()


# In[18]:


# If the same strain is present more than once, average its values
data = original_data.groupby('orf')[['chr5_ratio_wt_unt','chr1_ratio_wt_unt']].mean()


# In[19]:


data.drop(index='WT', inplace=True)


# In[20]:


print('Final data dimensions: %d x %d' % (data.shape))


# In[21]:


data.head()


# # Prepare the final dataset

# In[22]:


dataset_ids = [16447, 16446]


# In[23]:


datasets = datasets.reindex(index=dataset_ids)


# In[24]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[25]:


data.head()


# ## Subset to the genes currently in SGD

# In[26]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[27]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[28]:


data.head()


# # Normalize

# In[29]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[30]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[31]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[32]:


data.head()


# # Print out

# In[33]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[34]:


from IO.save_data_to_db3 import *


# In[35]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




