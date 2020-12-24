#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 31416893
paper_name = 'bhat_sengupta_2019' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/whole screen.xlsx', sheet_name='Sheet2')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data['Name'] = original_data['Name'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['Name'] = clean_genename(original_data['Name'])


# In[9]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['Name'], to='orf')


# In[10]:


name_fix_map = {'GON2': 'YLL033W','CRS5':'YOR031W','SRF5':'YOR041C','MOR1':'YDR366C','BOP1':'YPL221W','FMO': 'YHR176W','GON3':'YHR177W','FMP53':'YLR201C',
               'OCT':'YKL134C','FMP17':'YGR033C','FMP31':'YOR286W','SRF6':'YNL179C','SWS1':'YDR290W','AAD6':'YFL056C','YSN1':'YNR065C','FMP35':'YIL157C','SDL1':'YIL167W',
                'ZSP1':'YBR287W','SRF4':'YDL023C','FLO8':'YER109C'}


# In[11]:


for g in name_fix_map.keys():
    original_data.loc[original_data['orfs']==g,'orfs'] = name_fix_map[g]


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[13]:


original_data['data'] = original_data['Cysteine'] / original_data['Control']


# In[15]:


original_data.set_index('orfs', inplace=True)
original_data.index.name='orf'


# In[16]:


original_data = original_data.groupby(original_data.index).mean()


# In[17]:


original_data.shape


# # Prepare the final dataset

# In[18]:


dataset_ids = [16547]


# In[19]:


datasets = datasets.reindex(index=dataset_ids)


# In[22]:


data = original_data[['data']].copy()


# In[23]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[24]:


data.head()


# ## Subset to the genes currently in SGD

# In[27]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[28]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[29]:


data.head()


# In[31]:


data.shape


# # Normalize

# In[30]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[32]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[33]:


data_norm[data.isnull()] = np.nan


# In[34]:


data_all = data.join(data_norm)


# In[35]:


data_all.head()


# # Print out

# In[36]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[37]:


from IO.save_data_to_db3 import *


# In[38]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




