#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 26149385
paper_name = 'sauerwald_rapaport_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data (1)

# In[5]:


original_data1 = pd.read_excel('raw_data/zmb999100949sd2.xlsx', sheet_name='Supplementary Table 1')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))


# In[7]:


original_data1.head()


# In[8]:


original_data1['orf'] = original_data1['ORF']


# In[9]:


# Eliminate all white spaces & capitalize
original_data1['orf'] = clean_orf(original_data1['orf'])


# In[10]:


# Translate to ORFs 
original_data1['orf'] = translate_sc(original_data1['orf'], to='orf')


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orf'])
print(original_data1.loc[~t,])


# In[12]:


original_data1['data'] = -1


# In[13]:


original_data1.set_index('orf', inplace=True)


# In[14]:


original_data1 = original_data1[['data']].copy()


# In[15]:


original_data1 = original_data1.groupby(original_data1.index).mean()


# In[16]:


original_data1.shape


# # Load & process the data (2)

# In[17]:


original_data2 = pd.read_excel('raw_data/zmb999100949sd3.xlsx', sheet_name='Supplemental Table 2')


# In[18]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[19]:


original_data2.head()


# In[20]:


original_data2['orf'] = original_data2['ORF'].astype(str)


# In[21]:


original_data2['orf'] = clean_orf(original_data2['orf'])


# In[22]:


original_data2.loc[original_data2['orf'] == 'YOLO62C','orf'] = 'YOL062C'
original_data2.loc[original_data2['orf'] == 'YOLO57W','orf'] = 'YOL057W'


# In[23]:


# Translate to ORFs 
original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')


# In[24]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[25]:


original_data2['data'] = original_data2['Ratio SGal-Ura/ SGal complete'] - 1


# In[26]:


original_data2.set_index('orf', inplace=True)


# In[27]:


original_data2 = original_data2[['data']].copy()


# In[28]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[29]:


original_data2.shape


# # Merge

# In[45]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[46]:


original_data[original_data.isnull()] = 0


# In[47]:


original_data.head()


# # Prepare the final dataset

# In[48]:


data = original_data.copy()


# In[49]:


dataset_ids = [11865, 11866]
datasets = datasets.reindex(index=dataset_ids)


# In[50]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[51]:


data.head()


# ## Subset to the genes currently in SGD

# In[52]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[53]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[54]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[55]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[56]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[57]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[58]:


from IO.save_data_to_db3 import *


# In[59]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




