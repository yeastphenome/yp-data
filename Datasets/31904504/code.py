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


paper_pmid = 31904504
paper_name = 'zhao_deng_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/TableS2.xlsx', sheet_name='Sheet1', skiprows=2)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[8]:


original_data.columns = ['orf','gene','cfu','spot']


# In[9]:


original_data['orf'] = original_data['orf'].astype(str)


# In[10]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[11]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[13]:


original_data = original_data.loc[t,:]


# In[14]:


original_data['data'] = -1


# In[16]:


original_data.set_index('orf', inplace=True)


# # Prepare the final dataset

# In[18]:


dataset_ids = [16433]


# In[19]:


datasets = datasets.reindex(index=dataset_ids)


# In[20]:


data = original_data[['data']].copy()


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

# In[25]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[26]:


from IO.save_data_to_db2 import *


# In[27]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[28]:


save_data_to_db(data, paper_pmid)


# In[ ]:




