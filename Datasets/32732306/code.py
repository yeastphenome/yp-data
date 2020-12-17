#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd

import sys

from os.path import expanduser
sys.path.append(expanduser('~') + '/Lab/Utils/Python/')

from Conversions.translate import *
from Strings.is_a import *
from Math.normalize import z_transform_mode


# # Initial setup

# In[2]:


paper_pmid = 32732306
paper_name = 'silva_ideker_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# In[5]:


path_to_genes = '../../Private-Utils/datasets_gene2.txt'
path_to_consensus_tested = '../../Private-Utils/yp_2020-09-01_orfs.txt'


# # Load & process the data

# In[6]:


original_data = pd.read_csv('raw_data/File_S3_data.txt', header=None, names=['orfs','colony_unt','colony_trt','lagvstall_unt','lagvstall_trt'], sep='\t')


# In[7]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[8]:


original_data['orfs'] = original_data['orfs'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[10]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[12]:


original_data = original_data.loc[t,]


# In[13]:


original_data.set_index('orfs', inplace=True)


# In[14]:


original_data.index.name='orf'


# In[15]:


original_data['colony'] = original_data['colony_trt'] - original_data['colony_unt']
original_data['lagvstall'] = original_data['lagvstall_trt'] - original_data['lagvstall_unt']


# In[16]:


original_data = original_data.groupby(original_data.index).mean()


# In[17]:


original_data.shape


# # Prepare the final dataset

# In[18]:


data = original_data[['lagvstall','colony']].copy()


# In[19]:


dataset_ids = [16612,16613]
datasets = datasets.reindex(index=dataset_ids)


# In[20]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[21]:


print('Final data dimensions: %d x %d' % (data.shape))


# In[22]:


## Subset to the genes currently in SGD


# In[23]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=original_data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[24]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# # Normalize

# In[25]:


def normalize_phenotypic_scores(df, has_tested=False):
    
    if not has_tested:
        
        genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
        genes = genes.reset_index().set_index('systematic_name', drop=False)
        
        yp_orfs = pd.read_csv(path_to_consensus_tested, header=None)
        yp_orfs = yp_orfs[1].values
        consensus_tested = [tuple(genes.loc[orf,['id','systematic_name']]) for orf in yp_orfs]
        
        df_index = [tuple(x) for x in df.index]
        
        consensus_tested = list(set(consensus_tested + df_index))
        
        df = df.reindex(index=consensus_tested, fill_value=0)
        
    df_norm = z_transform_mode(df)
    
    return df_norm


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


# In[29]:


data_all = data.join(data_norm)


# In[30]:


data_all.head()


# # Print out

# In[31]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[32]:


from IO.save_data_to_db3 import *


# In[33]:


save_data_to_db(data_all, paper_pmid)

