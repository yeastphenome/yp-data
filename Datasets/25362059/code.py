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


paper_pmid = 25362059
paper_name = 'uehara_shimoi_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[8]:


original_data = pd.read_csv('raw_data/HEMF_production_capacity.txt', header=None, names=['orf','data'], sep='\t')


# In[9]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[11]:


original_data['orf'] = original_data['orf'].astype(str)


# In[12]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[13]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[14]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[15]:


original_data.set_index('orf', inplace=True)


# In[16]:


# These data represet a production capacity ratio (Hmut/Hwt). Shift data by 1 so that the WT phenotype is 0 (by convention)
original_data = original_data-1


# # Prepare the final dataset

# In[18]:


dataset_ids = [16495]


# In[19]:


datasets = datasets.reindex(index=dataset_ids)


# In[20]:


data = original_data.copy()


# In[21]:


data.columns = datasets['name'].values


# In[22]:


data = data.groupby(data.index).mean()


# In[23]:


# Create row index
data.index.name='orf'


# In[24]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[26]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[27]:


from IO.save_data_to_db2 import *


# In[28]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[29]:


save_data_to_db(data, paper_pmid)


# In[ ]:




