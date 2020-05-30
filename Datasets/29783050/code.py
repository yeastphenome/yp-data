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


paper_pmid = 29783050
paper_name = 'jiang_papadopoulos_2018' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[13]:


original_data = pd.read_csv('raw_data/hits.txt', sep='\t', header=None, names=['orfs'])


# In[14]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[15]:


original_data.head()


# In[16]:


original_data['orfs'] = original_data['orfs'].astype(str)


# In[17]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[18]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[19]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])


# In[20]:


print(original_data.loc[~t,])


# In[21]:


original_data['score'] = -1


# In[22]:


# If the same strain is present more than once, average its values
data = original_data.groupby('orfs')['score'].mean().to_frame()


# In[23]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Prepare the final dataset

# In[24]:


dataset_ids = [16444]


# In[25]:


datasets = datasets.reindex(index=dataset_ids)


# In[28]:


data.columns = datasets['name'].values


# In[29]:


# Create row index
data.index.name='orf'


# # Print out

# In[31]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[21]:


from IO.save_data_to_db2 import *


# In[ ]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[22]:


save_data_to_db(data, paper_pmid)

