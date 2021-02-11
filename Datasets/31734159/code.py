#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 31734159
paper_name = 'kuroda_avalos_2019' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/1-s2.0-S2405471219303825-mmc2.xlsx', sheet_name='1st screen', skiprows=4)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


hit_orfs = np.hstack([original_data.iloc[:,c].values for c in [0,5,10,19]])


# In[9]:


hit_orfs.shape


# In[10]:


hit_orfs = hit_orfs.astype(str)


# In[11]:


hit_data1 = np.hstack([original_data.iloc[:,c].values for c in [1,6,11,20]])
hit_data2 = np.hstack([original_data.iloc[:,c].values for c in [3,8,13]]+[np.zeros((original_data.shape[0]))+np.nan])


# In[12]:


hit_data1.shape


# In[13]:


hit_data2.shape


# In[14]:


original_data2 = pd.DataFrame(index=hit_orfs, columns=['0_iso','tf'], data=np.vstack((hit_data1, hit_data2)).T)


# In[15]:


original_data2['orfs'] = original_data2.index.values.astype(str)


# In[16]:


# Eliminate all white spaces & capitalize
original_data2['orfs'] = clean_orf(original_data2['orfs'])


# In[17]:


typo_fix = {'VER093C-A': 'YER093C-A','YMROB4W':'YMR084W','YARD02C-A':'YAR002C-A','YFLOOIW':'YFL001W',
           'YMLOIOC-B':'YML010C-B','VER091C':'YER091C','VCR086W':'YCR086W','YNLO15W':'YNL015W',
           'VER064C':'YER064C','VBR285W':'YBR285W','YMR08IC':'YMR081C','YJRIOOC':'YJR100C','YGR2B8W':'YGR288W',
           'YPROT4C':'YPR014C','VCR087W':'YCR087W','YARD50W':'YAR050W','VCL075W':'YCL075W','YJR09IC':'YJR091C',
           'YJR044C6':'YJR044C','VCR031C':'YCR031C'}


# In[18]:


original_data2['orfs'] = original_data2['orfs'].apply(lambda x: typo_fix[x] if x in list(typo_fix.keys()) else x)


# In[19]:


original_data2.shape


# In[20]:


original_data2 = original_data2.groupby('orfs').mean().reset_index()


# In[21]:


original_data2.shape


# In[22]:


original_data2.head()


# In[23]:


# Translate to ORFs 
original_data2['orfs'] = translate_sc(original_data2['orfs'], to='orf')


# In[24]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[25]:


original_data2 = original_data2.loc[t,:]


# In[26]:


original_data2.shape


# In[27]:


original_data2 = original_data2.groupby('orfs').mean()


# In[28]:


original_data2.head()


# In[29]:


original_data2.index.name='orf'


# In[30]:


original_data2.shape


# In[37]:


original_data2.loc['YBR010W',]


# # Prepare the final dataset

# In[38]:


data = original_data2[['0_iso','tf']].copy()


# In[39]:


dataset_ids = [16414,16411]
datasets = datasets.reindex(index=dataset_ids)


# In[40]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[41]:


data.head()


# ## Subset to the genes currently in SGD

# In[42]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[43]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[44]:


data.head()


# # Normalize

# In[45]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[46]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[47]:


data_norm[data.isnull()] = np.nan


# In[48]:


data_all = data.join(data_norm)


# In[49]:


data_all.head()


# # Print out

# In[50]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[51]:


from IO.save_data_to_db3 import *


# In[52]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




