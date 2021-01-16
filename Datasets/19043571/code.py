#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 19043571
paper_name = 'yu_bellaoui_2008' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[6]:


original_data = pd.read_excel('raw_data/13_15_data.xlsx', sheet_name='13&15 diploid (Figure 2)')


# In[7]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[8]:


original_data.head()


# In[9]:


original_data['orf'] = original_data['strain'].apply(lambda x: x.split(':')[0])


# In[10]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[11]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[13]:


# Split het and hom
original_data1 = original_data.loc[original_data['zygosity']=='hom'].copy()
original_data2 = original_data.loc[original_data['zygosity']=='het'].copy()


# In[14]:


original_data1.set_index('orf', inplace=True)
original_data2.set_index('orf', inplace=True)


# In[15]:


data_cols = ['z_result_nq:07_04_18_t01:cpmd13:0.98:ug/ml::::5_20:heho_05_06:YPD',
             'z_result_nq:07_04_18_t02:cmpd15:1:ug/ml::::5_20:heho_05_06:YPD']
original_data1 = original_data1[data_cols].copy()
original_data2 = original_data2[data_cols].copy()


# In[16]:


original_data1 = original_data1.apply(pd.to_numeric, axis=1, errors='coerce')
original_data2 = original_data2.apply(pd.to_numeric, axis=1, errors='coerce')


# In[17]:


original_data1 = original_data1.groupby(original_data1.index).mean()
original_data2 = original_data2.groupby(original_data2.index).mean()


# In[18]:


print(original_data1.shape)
print(original_data2.shape)


# # Load (haploid)

# In[19]:


original_data3 = pd.read_excel('raw_data/HAP compendium May 12 19 2006 including compounds 13 and 15-without badtags.xlsx', 
                               sheet_name='Sheet1')
print('Original data dimensions: %d x %d' % (original_data3.shape))


# In[21]:


original_data3.head()


# In[22]:


original_data3['orf'] = original_data3['strain'].astype(str)


# In[23]:


original_data3['orf'] = clean_orf(original_data3['orf'])


# In[24]:


original_data3['orf'] = translate_sc(original_data3['orf'], to='orf')


# In[25]:


t = looks_like_orf(original_data3['orf'])
print(original_data3.loc[~t,])


# In[26]:


original_data3.set_index('orf', inplace=True)
original_data3 = original_data3[['13Cmpd','15Cmpd']].astype(float)


# In[27]:


original_data3 = original_data3.groupby(original_data3.index).mean()


# In[28]:


print(original_data3.shape)


# # Merge

# In[34]:


original_data = pd.concat([original_data1,original_data2,original_data3], axis=1)


# In[36]:


# Flipping the sign on all because fitness scores are originally reported on the UNT/TRT scale
original_data = -original_data


# In[38]:


original_data.index.name='orf'


# # Prepare the final dataset

# In[40]:


data = original_data.copy()


# In[41]:


dataset_ids = [508,4982,5007,5008,511,4983]
datasets = datasets.reindex(index=dataset_ids)


# In[42]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[43]:


data.head()


# ## Subset to the genes currently in SGD

# In[44]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[45]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[46]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[47]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[48]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[49]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[50]:


from IO.save_data_to_db3 import *


# In[51]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




