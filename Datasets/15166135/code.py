#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 15166135
paper_name = 'lesage_bussey_2004' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/Table3-4.xlsx', sheet_name='Table 3', skiprows=2)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data['orf'] = original_data['ORF'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[10]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[12]:


original_data = original_data.loc[t,:]


# In[13]:


original_data.set_index('orf', inplace=True)


# In[14]:


# Separate hap from het data
original_data1 = original_data.loc[original_data['Haploid'].isnull(),]
original_data2 = original_data.loc[original_data['Haploid'].notnull(),]


# In[15]:


original_data1['data'] = original_data1['Heterozygous diploid'].apply(lambda x: -len(x.strip()))


# In[16]:


original_data2['data'] = original_data2['Haploid'].apply(lambda x: -len(x.strip()))


# In[17]:


original_data1 = original_data1[['data']].copy()
original_data2 = original_data2[['data']].copy()


# In[18]:


original_data1 = original_data1.groupby(original_data1.index).mean()
original_data2 = original_data2.groupby(original_data2.index).mean()


# # Load data (2)

# In[19]:


original_data3 = pd.read_excel('raw_data/Table3-4.xlsx', sheet_name='Table 4')


# In[20]:


original_data3.head()


# In[21]:


original_data3['orf'] = original_data3['ORF'].astype(str)


# In[22]:


# Eliminate all white spaces & capitalize
original_data3['orf'] = clean_orf(original_data3['orf'])


# In[23]:


# Translate to ORFs 
original_data3['orf'] = translate_sc(original_data3['orf'], to='orf')


# In[24]:


# Make sure everything translated ok
t = looks_like_orf(original_data3['orf'])
print(original_data3.loc[~t,])


# In[25]:


original_data3 = original_data3.loc[t,:]


# In[26]:


original_data3.set_index('orf', inplace=True)


# In[27]:


original_data3['data'] = 1


# In[28]:


original_data3 = original_data3[['data']].copy()


# In[29]:


original_data3 = original_data3.groupby(original_data3.index).mean()


# # Merge

# In[30]:


original_data = pd.concat([original_data2, original_data3], axis=0)


# In[33]:


original_data = pd.concat([original_data, original_data1], axis=1)


# In[37]:


original_data.shape


# In[44]:


original_data.index.name='orf'


# # Prepare the final dataset

# In[45]:


data = original_data.copy()


# In[46]:


dataset_ids = [21846,21848]
datasets = datasets.reindex(index=dataset_ids)


# In[47]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[48]:


data.head()


# ## Subset to the genes currently in SGD

# In[49]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[50]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[51]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[52]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[53]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[54]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[55]:


from IO.save_data_to_db3 import *


# In[56]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




