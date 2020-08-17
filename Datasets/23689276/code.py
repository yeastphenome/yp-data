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


paper_pmid = 23689276
paper_name = 'bowie_fyles_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[12]:


original_data1 = pd.read_excel('raw_data/c3ob40593a.xlsx', sheet_name='Table S1', skiprows=2)
original_data2 = pd.read_excel('raw_data/c3ob40593a.xlsx', sheet_name='Table S2', skiprows=2)


# In[9]:


print('Original data dimensions: %d x %d' % (original_data1.shape))


# In[13]:


original_data1['Orf'] = original_data1['Orf'].astype(str)
original_data2['Orf'] = original_data2['Orf'].astype(str)


# In[14]:


# Eliminate all white spaces & capitalize
original_data1['Orf'] = clean_orf(original_data1['Orf'])
original_data2['Orf'] = clean_orf(original_data2['Orf'])


# In[15]:


# Translate to ORFs 
original_data1['Orf'] = translate_sc(original_data1['Orf'], to='orf')
original_data2['Orf'] = translate_sc(original_data2['Orf'], to='orf')


# In[17]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['Orf'])
print(original_data1.loc[~t,])


# In[18]:


original_data1 = original_data1.loc[t,:]


# In[19]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['Orf'])
print(original_data2.loc[~t,])


# In[20]:


original_data2 = original_data2.loc[t,:]


# In[21]:


original_data1.set_index('Orf', inplace=True)
original_data2.set_index('Orf', inplace=True)


# In[24]:


data = original_data1['Ratio'].to_frame().join(original_data2['Ratio'].to_frame(), lsuffix='_1', rsuffix='_2')


# In[25]:


data.head()


# # Prepare the final dataset

# In[26]:


dataset_ids = [16557, 16558]


# In[27]:


datasets = datasets.reindex(index=dataset_ids)


# In[28]:


data.columns = datasets['name'].values


# In[33]:


for c in data.columns:
    data[c] = pd.to_numeric(data[c], errors='coerce')


# In[34]:


data = data.groupby(data.index).mean()


# In[35]:


# Create row index
data.index.name='orf'


# In[36]:


print('Final data dimensions: %d x %d' % (data.shape))


# In[37]:


data.head()


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




