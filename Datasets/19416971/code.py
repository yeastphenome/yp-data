#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 19416971
paper_name = 'khozoie_avery_2009' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/jbc.M109.005843-1.xls', sheet_name='Initial QN Screen', skiprows=3)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[10]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'], to='orf')


# In[11]:


to_drop = original_data.loc[original_data['ORF']=='EMPTY',]


# In[12]:


original_data.drop(index=to_drop.index, inplace=True)


# In[13]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[14]:


original_data = original_data.loc[t,:]


# In[15]:


# Eliminate the slow growin strains for which no accurate growth ratio could be calculated
to_drop = original_data.loc[original_data['Adjusted GR']=='SLOW',]


# In[16]:


original_data.drop(index=to_drop.index, inplace=True)


# In[17]:


original_data.shape


# In[18]:


original_data['Growth Ratio (GR)'] = pd.to_numeric(original_data['Growth Ratio (GR)'], errors='coerce')


# In[19]:


# Reverse the growth ratio so that lower values correspond to decreased growth and viceversa (originally, GR is reported as untreated vs treated)
original_data['GR2'] = 1 / original_data['Growth Ratio (GR)']


# In[20]:


# Normalize by plate median (as done oridinally)
def normalize_by_plate_median(plate_data):
    plate_median = plate_data['GR2'].median()
    plate_data['GR2_adjusted'] = plate_data['GR2'] / plate_median
    return plate_data


# In[21]:


original_data2 = original_data.groupby('Plate').apply(normalize_by_plate_median)


# In[22]:


original_data2.head()


# In[23]:


original_data2.set_index('ORF', inplace=True)
original_data2.index.name='orf'


# In[24]:


original_data2['data'] = original_data2['GR2_adjusted']


# In[25]:


original_data2 = original_data2[['data']].copy()


# In[26]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[27]:


original_data2.head()


# In[28]:


original_data2.shape


# # Prepare the final dataset

# In[29]:


data = original_data2.copy()


# In[30]:


dataset_ids = [16533]
datasets = datasets.reindex(index=dataset_ids)


# In[31]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[32]:


data.head()


# ## Subset to the genes currently in SGD

# In[33]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[34]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[35]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[36]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[37]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[38]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[39]:


from IO.save_data_to_db3 import *


# In[40]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




