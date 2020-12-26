#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 28520786
paper_name = 'choi_oshea_2017' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[12]:


xl = pd.ExcelFile('raw_data/journal.pone.0176085.s004.xlsx')
sheet_names = np.array(xl.sheet_names)[[0,2]]


# In[13]:


sheet_names


# In[18]:


original_data_list = []
for s in sheet_names:
    original_data = xl.parse(s)
    print('Original data dimensions: %d x %d' % (original_data.shape))
    original_data['orf'] = original_data['Gene_systematic name'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data.loc[original_data['orf']=='YLR287-A','orf'] = 'YLR287C-A'
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data['data'] = original_data.iloc[:,2].astype(float)
    original_data.set_index('orf', inplace=True)
    original_data = original_data[['data']].copy()
    original_data = original_data.groupby(original_data.index).mean()
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[19]:


original_data1, original_data2 = original_data_list


# In[20]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[25]:


# M&M indicates that the list includes 4974 knockout alleles of
# nonessential genes and 878 hypomorphic alleles of essential genes. Must
# eliminate the essential genes (no other way available)

noness = ~is_essential(original_data.index.values)
original_data = original_data.loc[noness.values,:]
print(original_data.shape)


# In[26]:


dataset_ids = [11807, 11808]


# # Prepare the final dataset

# In[27]:


data = original_data.copy()


# In[28]:


datasets = datasets.reindex(index=dataset_ids)


# In[29]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[30]:


data.head()


# ## Subset to the genes currently in SGD

# In[31]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[32]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[33]:


data.head()


# # Normalize

# In[34]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[35]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[36]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[37]:


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




