#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 18780730
paper_name = 'sinha_steinmetz_2008' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_csv('raw_data/deletion_pool_data.txt', sep='\t')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.columns


# In[8]:


original_data.head()


# In[9]:


original_data['orfs'] = original_data['orf::batch:tagtype'].apply(lambda x: x.split(':')[0])


# In[10]:


original_data['orfs'] = original_data['orfs'].astype(str)


# In[11]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[12]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[13]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[14]:


original_data['37C'] = (original_data['37C_T5'] / original_data['T0']) / (original_data['30C_T5'] / original_data['T0'])


# In[15]:


original_data.sort_values(by='37C', ascending=False)[['orfs','T0','30C_T1','30C_T2','30C_T3','30C_T4','30C_T5','37C_T1','37C_T2','37C_T3','37C_T4','37C_T5']].head()


# In[16]:


original_data['rapa_12h'] = original_data['30C_RAPA_T1'] / original_data['T0']
original_data['rapa_24h'] = original_data['30C_RAPA_T2'] / original_data['T0']
original_data['rapa_36h'] = original_data['30C_RAPA_T3'] / original_data['T0']


# In[17]:


original_data.set_index('orfs', inplace=True)
original_data.index.name='orf'


# In[18]:


# Splits homozygous and heterozygous mutants
original_data_hom = original_data.loc[original_data['zygosity']=='hom'].copy()
original_data_het = original_data.loc[original_data['zygosity']=='het'].copy()


# In[19]:


original_data_hom = original_data_hom.groupby(original_data_hom.index).mean()
original_data_hom.shape


# In[20]:


original_data_het = original_data_het.groupby(original_data_het.index).mean()
original_data_het.shape


# In[21]:


original_data_hom = original_data_hom[['37C','rapa_12h','rapa_24h','rapa_36h']].copy()
original_data_het = original_data_het[['37C','rapa_12h','rapa_24h','rapa_36h']].copy()


# In[22]:


original_data_hom = original_data_hom.groupby(original_data_hom.index).mean()
original_data_het = original_data_het.groupby(original_data_het.index).mean()


# In[23]:


# Pull them back together
original_data = original_data_hom.join(original_data_het, lsuffix='_hom', rsuffix='_het', how='outer')


# In[24]:


original_data.shape


# In[25]:


original_data.head()


# # Prepare the final dataset

# In[26]:


data = original_data.copy()


# In[27]:


dataset_ids = [16511,16639,16640,16641,16638,16642,16643,16644]
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




