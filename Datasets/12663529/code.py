#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')

import re


# # Initial setup

# In[2]:


paper_pmid = 12663529
paper_name = 'page_bussey_2003' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


files = ['Table 2.txt','Table 3.txt','Table 4.txt']
patterns = ['(wt|\d+)\s(wt|\d+)\s(wt|\d+|NA)','(wt|\s\d+)\s(wt|\d+)\s(wt|\d+|NA)','\s(wt|\d+)']


# In[6]:


original_data_list = []
for ixf, f in enumerate(files):
    original_data = pd.read_csv('raw_data/' + f, header=None, sep='\n')
    print('Original data dimensions: %d x %d' % (original_data.shape))
    print(original_data.head())
    original_data['orf'] = original_data[0].apply(lambda x: x.split(' ')[1] if looks_like_orf(x.split(' ')[1]) else x.split(' ')[0])
    original_data['info'] = original_data[0].apply(lambda x: re.search(patterns[ixf], x))
    original_data['info2'] = original_data['info'].apply(lambda x: x.group(0) if x else 'nan nan nan')
    original_data['info2'] = original_data['info2'].apply(lambda x: [y for y in x.strip().split(' ')])
    original_data[['data1','data2','data3']] = pd.DataFrame(original_data['info2'].to_list())
    
    for c in ['data1','data2','data3']:
        original_data.loc[original_data[c]=='wt', c] = '100'
        original_data[c] = pd.to_numeric(original_data[c], errors='coerce')
        
    original_data['orf'] = original_data['orf'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data = original_data.loc[t,:]
    original_data.set_index('orf', inplace=True)
    original_data = original_data[['data1','data2','data3']].copy()
    original_data = original_data.groupby(original_data.index).mean()
    
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[7]:


original_data = pd.concat(original_data_list, axis=0)


# In[8]:


original_data.shape


# In[9]:


original_data = 100 - original_data


# In[10]:


original_data = original_data.groupby(original_data.index).mean()


# In[11]:


original_data.shape


# # Prepare the final dataset

# In[12]:


data = original_data.copy()


# In[13]:


dataset_ids = [81, 83, 82]
datasets = datasets.reindex(index=dataset_ids)


# In[14]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[15]:


data.head()


# ## Subset to the genes currently in SGD

# In[16]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[17]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[18]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[19]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[20]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[21]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[22]:


from IO.save_data_to_db3 import *


# In[23]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




