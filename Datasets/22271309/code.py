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


paper_pmid = 22271309
paper_name = 'yibmantasiri_bellows_2012' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, names=['genes'])


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data['genes'] = original_data['genes'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[9]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[11]:


original_data['data'] = -1


# In[12]:


original_data.set_index('orfs', inplace=True)


# In[18]:


original_data.head()


# # Prepare the final dataset

# In[13]:


dataset_ids = [16485]


# In[14]:


datasets = datasets.reindex(index=dataset_ids)


# In[26]:


data = original_data['data'].copy()


# In[27]:


data.columns = datasets['name'].values


# In[28]:


data = data.groupby(data.index).mean().to_frame()


# In[29]:


# Create row index
data.index.name='orf'


# In[31]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[32]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[33]:


from IO.save_data_to_db2 import *


# In[34]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[35]:


save_data_to_db(data, paper_pmid)


# In[ ]:



