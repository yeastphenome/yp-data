#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 21685046
paper_name = 'yu_vitek_2011' 


# In[39]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[40]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data (hap)

# In[5]:


original_data = pd.read_csv('raw_data/PiiMS_Dataset_Full_genome_knockout_haploid_1469056539143/Dataset_387_Aggregate_1299870873410.csv')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data['orf'] = original_data['Line'].apply(lambda x: x.split('_')[0])


# In[9]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[12]:


original_data.loc[original_data['orf']=='YLR287-A','orf'] = 'YLR287C-A'


# In[13]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[14]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[15]:


original_data.set_index('orf', inplace=True)


# In[16]:


original_data.drop(columns=['Line'], inplace=True)


# In[18]:


for c in original_data.columns:
    original_data[c] = pd.to_numeric(original_data[c], errors='coerce')


# In[19]:


original_data = original_data.groupby(original_data.index).mean()


# In[20]:


original_data.shape


# # Load & process the data (het)

# In[21]:


original_data2 = pd.read_csv('raw_data/PiiMS_Dataset_Full_genome_knockout_diploid_1469056540949/Dataset_390_Aggregate_1299877732527.csv')


# In[22]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[24]:


original_data2.head()


# In[25]:


original_data2['orf'] = original_data2['Line'].apply(lambda x: x.split('_')[0])


# In[26]:


# Eliminate all white spaces & capitalize
original_data2['orf'] = clean_orf(original_data2['orf'])


# In[27]:


# Translate to ORFs 
original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')


# In[28]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[29]:


original_data2.set_index('orf', inplace=True)


# In[30]:


original_data2.drop(columns=['Line'], inplace=True)


# In[31]:


for c in original_data2.columns:
    original_data2[c] = pd.to_numeric(original_data2[c], errors='coerce')


# In[32]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[33]:


original_data2.shape


# # Merge

# In[34]:


original_data = original_data.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[35]:


original_data.head()


# # Prepare the final dataset

# In[42]:


data = original_data.copy()


# In[43]:


dataset_ids = np.arange(4797, 4825)
datasets = datasets.reindex(index=dataset_ids)


# In[44]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[45]:


data.head()


# ## Subset to the genes currently in SGD

# In[46]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[47]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[48]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[49]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[50]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[51]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[52]:


from IO.save_data_to_db3 import *


# In[53]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




