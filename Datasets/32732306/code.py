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


paper_pmid = 32732306
paper_name = 'silva_ideker_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[12]:


original_data = pd.read_csv('raw_data/File_S3_data.txt', header=None, names=['orfs','colony_unt','colony_trt','lagvstall_unt','lagvstall_trt'], sep='\t')


# In[13]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[14]:


original_data['orfs'] = original_data['orfs'].astype(str)


# In[15]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[16]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[17]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[18]:


original_data = original_data.loc[t,]


# In[19]:


original_data.set_index('orfs', inplace=True)


# In[24]:


original_data['colony'] = original_data['colony_trt'] - original_data['colony_unt']
original_data['lagvstall'] = original_data['lagvstall_trt'] - original_data['lagvstall_unt']


# # Prepare the final dataset

# In[27]:


data = original_data[['lagvstall','colony']].copy()


# In[29]:


dataset_ids = [16612,16613]


# In[30]:


datasets = datasets.reindex(index=dataset_ids)


# In[31]:


data.columns = datasets['name'].values


# In[32]:


data = data.groupby(data.index).mean()


# In[33]:


# Create row index
data.index.name='orf'


# In[34]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[37]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[38]:


from IO.save_data_to_db2 import *


# In[39]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[40]:


save_data_to_db(data, paper_pmid)


# In[ ]:




