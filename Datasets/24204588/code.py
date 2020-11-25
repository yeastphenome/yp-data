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


paper_pmid = 24204588
paper_name = 'darvishi_smith_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[15]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, names=['orfs','score'], sep='\t')


# In[16]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[17]:


original_data['orfs'] = original_data['orfs'].astype(str)


# In[18]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[19]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[20]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[21]:


original_data['data'] = -original_data['score']/100


# In[22]:


original_data.set_index('orfs', inplace=True)


# # Prepare the final dataset

# In[24]:


dataset_ids = [16514]


# In[25]:


datasets = datasets.reindex(index=dataset_ids)


# In[26]:


data = original_data[['data']].copy()


# In[27]:


data.columns = datasets['name'].values


# In[28]:


data = data.groupby(data.index).mean()


# In[29]:


# Create row index
data.index.name='orf'


# In[30]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[31]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[32]:


from IO.save_data_to_db2 import *


# In[33]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[34]:


save_data_to_db(data, paper_pmid)


# In[ ]:




