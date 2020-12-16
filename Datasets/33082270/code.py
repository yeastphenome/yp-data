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


paper_pmid = 33082270
paper_name = 'berg_brandl_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# In[5]:


path_to_genes = '../../Private-Utils/datasets_gene2.txt'
path_to_consensus_tested = '../../Private-Utils/yp_2020-09-01_orfs.txt'


# # Load & process the data

# In[6]:


original_data = pd.read_excel('raw_data/SupplementalFile2_AZCInitialScreenZscores.xlsx', 
                              sheet_name='SupplementalFile1_Zscores', skiprows=5)


# In[7]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[8]:


original_data['orfs'] = original_data['Systematic Name'].astype(str)


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


original_data['data'] = original_data['Z-Score']


# In[13]:


ix = is_essential(original_data['orfs'])


# In[14]:


original_data = original_data.loc[~ix.values]


# In[15]:


original_data.set_index('orfs', inplace=True)
original_data.index.name='orf'


# In[16]:


original_data = original_data.groupby(original_data.index).mean()


# In[17]:


original_data.shape


# # Prepare the final dataset

# In[18]:


data = original_data[['data']].copy()


# In[19]:


dataset_ids = [16665]
datasets = datasets.reindex(index=dataset_ids)


# In[20]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[21]:


data.head()


# ## Subset to the genes currently in SGD

# In[22]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=original_data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[23]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[24]:


data.head()


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
        
        # Merge consensus tested with current list of ORFs to make sure 
        # we don't miss any gene just because it's not in the current consensus
        consensus_tested = list(set(consensus_tested + df_index))
        
        df = df.reindex(index=consensus_tested, fill_value=0)
        
    df_norm = z_transform_mode(data)
    
    return df_norm


# In[26]:


data_norm = normalize_phenotypic_scores(data)


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
    df = data_all.xs('value', level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[32]:


from IO.save_data_to_db3 import *


# In[36]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




