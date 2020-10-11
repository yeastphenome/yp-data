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


# # Initial setup

# In[2]:


paper_pmid = 24360837
paper_name = 'hoepfner_movva_2014' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data - Benomyl

# In[5]:


original_data1 = pd.read_csv('large_files/raw_data/HOP_scores-benomyl.txt', sep='\t')
original_data2 = pd.read_csv('large_files/raw_data/HIP_scores-benomyl.txt', sep='\t')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[7]:


# Keep the sensitivity scores, not z-scores (z-score normalize each strain to its phenotype to all other compounds in the dataset)


# In[8]:


cols1 = [c for c in original_data1.columns.values if 'z-score' not in c]
cols2 = [c for c in original_data2.columns.values if 'z-score' not in c]


# In[9]:


original_data1 = original_data1.loc[:, cols1]
original_data2 = original_data2.loc[:, cols2]


# In[10]:


orf_col = 'Systematic Name'


# In[11]:


original_data1[orf_col] = original_data1[orf_col].astype(str)
original_data2[orf_col] = original_data2[orf_col].astype(str)


# In[12]:


# Eliminate all white spaces & capitalize
original_data1[orf_col] = clean_orf(original_data1[orf_col])
original_data2[orf_col] = clean_orf(original_data2[orf_col])


# In[13]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1[orf_col], to='orf')
original_data2['orfs'] = translate_sc(original_data2[orf_col], to='orf')


# In[14]:


original_data1.loc[original_data1['orfs'] == 'YBR160WAS','orfs'] = 'YBR160W'
original_data2.loc[original_data2['orfs'] == 'YBR160WAS','orfs'] = 'YBR160W'


# In[15]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[16]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[17]:


original_data1 = original_data1.loc[t,:]
original_data2 = original_data2.loc[t,:]


# In[18]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)


# In[19]:


original_data1['data'] = original_data1.mean(axis=1)
original_data2['data'] = original_data2.mean(axis=1)


# In[20]:


original_data = original_data1[['data']].join(original_data2[['data']], how='outer', rsuffix='_hop', lsuffix='_hip')


# In[21]:


dataset_ids = [1087, 16622]


# In[22]:


datasets = datasets.reindex(index=dataset_ids)


# In[23]:


data = original_data[['data_hop','data_hip']].copy()


# In[24]:


data.columns = datasets['name'].values


# In[25]:


data = data.groupby(data.index).mean()


# In[26]:


# Create row index
data.index.name='orf'


# In[27]:


print('Final data dimensions: %d x %d' % (data.shape))


# In[28]:


data_benomyl = data.copy()


# In[29]:


data_benomyl.to_csv(paper_name + '_benomyl.txt', sep='\t')


# # Save to DB

# In[30]:


from IO.save_data_to_db2 import *


# In[31]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[32]:


data.head()


# In[33]:


save_data_to_db(data, paper_pmid, delete=False)


# In[ ]:




