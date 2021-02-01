#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 12615994
paper_name = 'zewaii_huang_2003' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


files = ['0118Table3.xlsx','0118Table4.xlsx','0118Table5.xlsx']


# In[20]:


original_data1_list = []
original_data2_list = []
for ixf, f in enumerate(files):
    sr = 3 if ixf == 0 else 2
    original_data = pd.read_excel('raw_data/' + f, sheet_name='Sheet1', skiprows=sr)
    print('Original data dimensions: %d x %d' % (original_data.shape))
    print(original_data.head())
    
    original_data.columns = [c.strip() for c in original_data.columns]
    
    original_data['orf'] = original_data['ORF'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    
    original_data.loc[original_data['orf']=='YLR287-A','orf'] = 'YLR287C-A'
    
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data = original_data.loc[t,:]
    
    original_data['data'] = original_data['Score'].apply(lambda x: len(x.strip()) if '+' in x else -len(x.strip()))
    
    original_data.set_index('orf', inplace=True)
    
    original_data1 = original_data.loc[original_data['Cell type']=='a'].copy()
    original_data2 = original_data.loc[original_data['Cell type']=='a/a'].copy()
    
    original_data1 = original_data1[['data']].copy()
    original_data1 = original_data1.groupby(original_data1.index).mean()
    
    original_data2 = original_data2[['data']].copy()
    original_data2 = original_data2.groupby(original_data2.index).mean()
    
    original_data1_list.append(original_data1)
    original_data2_list.append(original_data2)


# In[21]:


original_data1 = pd.concat(original_data1_list, axis=0)
original_data2 = pd.concat(original_data2_list, axis=0)


# In[22]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[23]:


original_data.head()


# # Prepare the final dataset

# In[24]:


data = original_data.copy()


# In[25]:


dataset_ids = [73, 414]
datasets = datasets.reindex(index=dataset_ids)


# In[26]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[27]:


data.head()


# ## Subset to the genes currently in SGD

# In[28]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[29]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[30]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[31]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[32]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[33]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[34]:


from IO.save_data_to_db3 import *


# In[35]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




