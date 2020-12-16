#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd

import sys
import re

from os.path import expanduser
sys.path.append(expanduser('~') + '/Lab/Utils/Python/')

from Conversions.translate import *
from Strings.is_a import *
from Math.normalize import z_transform_mode


# # Initial setup

# In[2]:


paper_pmid = 16365294
paper_name = 'ohya_morishita_2005' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# In[5]:


path_to_genes = '../../Private-Utils/datasets_gene2.txt'
path_to_consensus_tested = '../../Private-Utils/yp_2020-09-01_orfs.txt'


# # Load & process the data

# In[6]:


original_data = pd.read_csv('raw_data/mutant_analysis_2011_10_20.tab', sep='\t')


# In[7]:


original_data.head()


# In[8]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[9]:


original_data['name'] = original_data['name'].astype(str)


# In[10]:


# Eliminate all white spaces & capitalize
original_data['name'] = clean_orf(original_data['name'])


# In[11]:


# Translate to ORFs 
original_data['name'] = translate_sc(original_data['name'], to='orf')


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data['name'])
print(original_data.loc[~t,])


# In[13]:


original_data.set_index('name', inplace=True)
original_data.index.name='orf'


# In[14]:


original_data = original_data.groupby(original_data.index).mean()


# In[15]:


original_data.shape


# # Process parameters

# In[16]:


datasets.head()


# In[17]:


datasets['phenotype_name'] = ''
datasets['phenotype_id'] = ''
for d in datasets.index.values:
    t = datasets.loc[d,'name'].split('|')[1]
    matches = re.findall('\(([A-Z0-9_\-]*?)\)', t)
    
    datasets.loc[d,'phenotype_name'] = t
    datasets.loc[d,'phenotype_id'] = matches[0]


# In[18]:


# Exclude the CV parameters and the ones that are not easily interpretable
original_data = original_data.reindex(columns=datasets['phenotype_id'].values)


# # Prepare the final dataset

# In[19]:


data = original_data.copy()


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
# 

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


# In[ ]:




