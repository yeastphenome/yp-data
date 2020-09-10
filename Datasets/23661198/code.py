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


paper_pmid = 23661198
paper_name = 'vahey_voldman_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[14]:


original_data1 = pd.read_excel('raw_data/c3lc50162k.xlsx', sheet_name='Sheet1', skiprows=1)
original_data2 = pd.read_excel('raw_data/c3lc50162k_2.xlsx', sheet_name='Sheet1', skiprows=1)


# In[15]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[17]:


original_data2.head()


# In[18]:


original_data1['Unnamed: 0'] = original_data1['Unnamed: 0'].astype(str)
original_data2['Unnamed: 0'] = original_data2['Unnamed: 0'].astype(str)


# In[19]:


# Eliminate all white spaces & capitalize
original_data1['Unnamed: 0'] = clean_orf(original_data1['Unnamed: 0'])
original_data2['Unnamed: 0'] = clean_orf(original_data2['Unnamed: 0'])


# In[20]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['Unnamed: 0'], to='orf')
original_data2['orfs'] = translate_sc(original_data2['Unnamed: 0'], to='orf')


# In[21]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[22]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[23]:


original_data1['data'] = original_data1['y (avg.)']
original_data2['data'] = original_data2['y (avg.)']


# In[24]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)


# In[29]:


data = original_data1['data'].to_frame().join(original_data2['data'].to_frame(), lsuffix='_1', rsuffix='_2', how='outer')


# In[30]:


data.head()


# # Prepare the final dataset

# In[31]:


dataset_ids = [95, 16597]


# In[32]:


datasets = datasets.reindex(index=dataset_ids)


# In[33]:


data.columns = datasets['name'].values


# In[34]:


data = data.groupby(data.index).mean()


# In[35]:


# Create row index
data.index.name='orf'


# In[36]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[38]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[39]:


from IO.save_data_to_db2 import *


# In[40]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[41]:


save_data_to_db(data, paper_pmid)


# In[ ]:




