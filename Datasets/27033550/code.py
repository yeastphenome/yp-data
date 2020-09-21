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


paper_pmid = 27033550
paper_name = 'yang_nystrom_2016' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[8]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, names=['orfs','values'], sep='\t')


# In[9]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[11]:


original_data['orfs'] = original_data['orfs'].astype(str)


# In[12]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[13]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[14]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[15]:


original_data['values'] = original_data['values'].astype(float)


# In[17]:


# Log-transform fold enrichments so that the mode (all non-hit strains) will be 0, by default
original_data['values'] = np.log2(original_data['values'])


# In[16]:


original_data.set_index('orfs', inplace=True)


# # Prepare the final dataset

# In[19]:


dataset_ids = [16413]


# In[20]:


datasets = datasets.reindex(index=dataset_ids)


# In[21]:


data = original_data.copy()


# In[22]:


data.columns = datasets['name'].values


# In[23]:


data = data.groupby(data.index).mean()


# In[24]:


# Create row index
data.index.name='orf'


# In[25]:


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




