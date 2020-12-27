#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')

from functools import reduce


# # Initial setup

# In[2]:


paper_pmid = 27188886
paper_name = 'pautasso_rossi_2016' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/Supplementary_Table_1.xlsx', sheet_name='ORFs 2<Z<2', skiprows=2, header=[0,1])


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data_list = []
a = ['TPK1','TPK2','TPK3','BCY1']
b = ['Activators','Inhibitors']
for ia in a:
    genes_list = []
    for ib in b:
        genes = original_data[(ia,ib)].to_frame()
        genes.columns = ['genes']
        genes['genes'] = genes['genes'].astype(str)
        genes['genes'] = clean_genename(genes['genes'])
        genes['orfs'] = translate_sc(genes['genes'], to='orf')
        t = looks_like_orf(genes['orfs'])
#         print(genes.loc[~t])
        genes = genes.loc[t,:]
        if ib == 'Activators':
            genes['data'] = -1
        else:
            genes['data'] = 1
        
        genes.set_index('orfs', inplace=True)
        genes.index.name='orf'
        genes = genes[['data']].copy()
        
        genes_list.append(genes)
            
    data = pd.concat(genes_list, axis=0, ignore_index=False)
    data = data.groupby(data.index).mean()
    
    original_data_list.append(data)


# In[18]:


original_data = reduce(lambda x, y: pd.merge(x, y, how='outer', on='orf'), original_data_list)


# In[20]:


original_data[original_data.isnull()] = 0


# In[21]:


original_data.head()


# In[22]:


original_data.shape


# # Prepare the final dataset

# In[23]:


data = original_data.copy()


# In[24]:


dataset_ids = [16247, 16248, 16249, 16250]
datasets = datasets.reindex(index=dataset_ids)


# In[25]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[26]:


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


# In[33]:


data_all.head()


# # Print out

# In[34]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[35]:


from IO.save_data_to_db3 import *


# In[36]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




