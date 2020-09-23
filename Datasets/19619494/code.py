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


paper_pmid = 19619494
paper_name = 'okamoto_ohsumi_2009' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# ### Part1

# In[13]:


original_data1 = pd.read_csv('raw_data/Table1.txt', header=None, names=['genes','data'], sep='\t')


# In[14]:


print('Original data dimensions: %d x %d' % (original_data1.shape))


# In[15]:


original_data1['genes'] = original_data1['genes'].astype(str)


# In[16]:


# Eliminate all white spaces & capitalize
original_data1['genes'] = clean_genename(original_data1['genes'])


# In[17]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['genes'], to='orf')


# In[18]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[20]:


data_dict = {'+++': 0, '++': -1, '+': -2, '-': -3}


# In[21]:


original_data1['data_score'] = original_data1['data'].apply(lambda x: data_dict[x])


# In[23]:


original_data1['data_score'] = original_data1['data_score'].astype(int)


# In[24]:


original_data1.set_index('orfs', inplace=True)


# ### Part2

# In[25]:


original_data2 = pd.read_csv('raw_data/TableS1.txt', header=None, names=['genes','data'], sep='\t')


# In[26]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[27]:


original_data2['genes'] = original_data2['genes'].astype(str)


# In[28]:


# Eliminate all white spaces & capitalize
original_data2['genes'] = clean_genename(original_data2['genes'])


# In[29]:


# Translate to ORFs 
original_data2['orfs'] = translate_sc(original_data2['genes'], to='orf')


# In[30]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[31]:


data_dict = {'+++': 0, '++': -1, '+': -2, '-': -3}


# In[32]:


original_data2['data_score'] = original_data2['data'].apply(lambda x: data_dict[x])


# In[33]:


original_data2['data_score'] = original_data2['data_score'].astype(int)


# In[34]:


original_data2.set_index('orfs', inplace=True)


# In[35]:


original_data2.head()


# In[38]:


data = original_data1['data_score'].to_frame().join(original_data2['data_score'].to_frame(), how='outer', lsuffix='_1', rsuffix='_s1')


# In[42]:


data['data_score'] = data.mean(axis=1)


# # Prepare the final dataset

# In[54]:


dataset_ids = [16606]


# In[55]:


datasets = datasets.reindex(index=dataset_ids)


# In[56]:


data = data['data_score'].to_frame()


# In[57]:


data.columns = datasets['name'].values


# In[58]:


data = data.groupby(data.index).mean()


# In[59]:


# Create row index
data.index.name='orf'


# In[60]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[61]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[62]:


from IO.save_data_to_db2 import *


# In[63]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[64]:


save_data_to_db(data, paper_pmid)


# In[ ]:




