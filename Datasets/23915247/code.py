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


paper_pmid = 23915247
paper_name = 'porcu_ragnini_wilson_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/12943_2012_1173_MOESM1_ESM.xls', sheet_name='Sheet1')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data['orfs'] = original_data['Feature Name'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[9]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[11]:


# Taking the opposite since the paper refers to sensitive hits as those with log2 > 0.5
original_data['data'] = -original_data['Average Log2 Ratio']


# In[12]:


original_data.set_index('orfs', inplace=True)


# # Prepare the final dataset

# In[13]:


dataset_ids = [16561]


# In[14]:


datasets = datasets.reindex(index=dataset_ids)


# In[15]:


data = original_data[['data']].copy()


# In[16]:


data.columns = datasets['name'].values


# In[17]:


data = data.groupby(data.index).mean()


# In[18]:


# Create row index
data.index.name='orf'


# In[19]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[20]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[21]:


from IO.save_data_to_db2 import *


# In[22]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[23]:


save_data_to_db(data, paper_pmid)


# In[ ]:



