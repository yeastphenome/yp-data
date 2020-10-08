#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import pandas as pd

import sys

from os.path import expanduser
sys.path.append(expanduser('~') + '/Lab/Utils/Python/')

from Conversions.translate import *
from Strings.is_a import *


# # Initial setup

# In[4]:


paper_pmid = 24035500
paper_name = 'shimada_gasser_2013' 


# In[5]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[6]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[7]:


original_data1 = pd.read_excel('raw_data/Shimada_etal_BHS345_HIPHOP_Data.xlsx', sheet_name='HIP')
original_data2 = pd.read_excel('raw_data/Shimada_etal_BHS345_HIPHOP_Data.xlsx', sheet_name='HOP')


# In[8]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[9]:


orf_col = 'SYSTEMATIC_NAME'


# In[10]:


original_data1[orf_col] = original_data1[orf_col].astype(str)
original_data2[orf_col] = original_data2[orf_col].astype(str)


# In[11]:


# Eliminate all white spaces & capitalize
original_data1[orf_col] = clean_orf(original_data1[orf_col])
original_data2[orf_col] = clean_orf(original_data2[orf_col])


# In[12]:


# Translate to ORFs 
original_data1[orf_col] = translate_sc(original_data1[orf_col], to='orf')
original_data2[orf_col] = translate_sc(original_data2[orf_col], to='orf')


# In[13]:


# Make sure everything translated ok
t = looks_like_orf(original_data1[orf_col])
print(original_data1.loc[~t,])


# In[14]:


# Make sure everything translated ok
t = looks_like_orf(original_data2[orf_col])
print(original_data2.loc[~t,])


# In[15]:


original_data1['data'] = original_data1['Z_SCORE']
original_data2['data'] = original_data2['Z_SCORE']


# In[16]:


original_data1.set_index(orf_col, inplace=True)
original_data2.set_index(orf_col, inplace=True)


# In[17]:


original_data = original_data1[['data']].join(original_data2[['data']], how='outer', lsuffix='_hip', rsuffix='_hop')


# # Prepare the final dataset

# In[20]:


dataset_ids = [16608, 16609]


# In[21]:


datasets = datasets.reindex(index=dataset_ids)


# In[23]:


data = original_data[['data_hip','data_hop']].copy()


# In[24]:


data.columns = datasets['name'].values


# In[25]:


data = data.groupby(data.index).mean()


# In[26]:


# Create row index
data.index.name='orf'


# In[27]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[30]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[31]:


from IO.save_data_to_db2 import *


# In[32]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[33]:


save_data_to_db(data, paper_pmid)


# In[ ]:




