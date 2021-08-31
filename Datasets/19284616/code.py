#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 19284616
paper_name = 'thorsen_tamas_2009' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


sheets1 = ['0.5 As 24','1 As 24h','1.5 As 24h','0.5As 48h','1 As 48h','1.5 As 48h']
sheets2 = ['Cd 75 24h','Cd 100 24h','Cd 150 24h','Cd 75 48h','Cd 100 48h','Cd 150 48h']


# In[6]:


original_data_list1 = []
for s in sheets1:
    original_data = pd.read_excel('raw_data/As(III) list avg.xlsx', sheet_name=s)
    print('Original data dimensions: %d x %d' % (original_data.shape))
#     print(original_data.head())
    original_data['orf'] = original_data['Sys name'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data = original_data.loc[t,:]
    original_data['data'] = original_data['Average M'].astype(float)
    original_data.set_index('orf', inplace=True)
    original_data = original_data[['data']].copy()
    original_data = original_data.groupby(original_data.index).mean()
    print(original_data.shape)
    original_data_list1.append(original_data)


# In[7]:


original_data_list2 = []
for s in sheets2:
    original_data = pd.read_excel('raw_data/Cd list avg.xlsx', sheet_name=s)
    print('Original data dimensions: %d x %d' % (original_data.shape))
#     print(original_data.head())
    original_data['orf'] = original_data['Sys name'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data = original_data.loc[t,:]
    original_data['data'] = original_data['Average M'].astype(float)
    original_data.set_index('orf', inplace=True)
    original_data = original_data[['data']].copy()
    original_data = original_data.groupby(original_data.index).mean()
    print(original_data.shape)
    original_data_list2.append(original_data)


# In[8]:


original_data1 = pd.concat(original_data_list1, axis=1)
original_data2 = pd.concat(original_data_list2, axis=1)


# In[9]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[16]:


# Despite the description in the Methods ("resistance and sensitivity is indicated by positive and negative values respectively"), 
# comparison to other published studies of Arsenite and Cadmium suggests that the opposite of the reported values should be taken.
original_data = -original_data


# In[17]:


original_data.shape


# In[18]:


original_data.head()


# # Prepare the final dataset

# In[20]:


data = original_data.copy()


# In[21]:


dataset_ids = [1320, 5362, 5363, 5364, 5365, 5366, 1319, 5367, 5368, 5369, 5370, 5371]
datasets = datasets.reindex(index=dataset_ids)


# In[22]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[23]:


data.head()


# ## Subset to the genes currently in SGD

# In[24]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[25]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[26]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[27]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[28]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[29]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[30]:


from IO.save_data_to_db3 import *


# In[31]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




