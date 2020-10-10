#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import pandas as pd

import sys

from os.path import expanduser
sys.path.append(expanduser('~') + '/Lab/Utils/Python/')

from Conversions.translate import *
from Strings.is_a import *


# # Initial setup

# In[4]:


paper_pmid = 24034557
paper_name = 'vandenbosch_coenye_2013' 


# In[5]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[6]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[7]:


original_data = pd.read_excel('raw_data/fyr12071-sup-0002-TableS1.xlsx', sheet_name='Blad1', skiprows=1)


# In[8]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[9]:


original_data.head()


# In[10]:


original_data['orfs'] = original_data['Unnamed: 0'].astype(str)


# In[11]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[12]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[13]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[14]:


original_data.set_index('orfs', inplace=True)


# # Prepare the final dataset

# In[15]:


dataset_ids = [16620, 16621]


# In[16]:


datasets = datasets.reindex(index=dataset_ids)


# In[17]:


data = original_data[['value','value.1']].copy()


# In[18]:


data.columns = datasets['name'].values


# In[19]:


data = data.groupby(data.index).mean()


# In[20]:


# Create row index
data.index.name='orf'


# In[21]:


print('Final data dimensions: %d x %d' % (data.shape))


# In[22]:


data.head()


# # Print out

# In[23]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[24]:


from IO.save_data_to_db2 import *


# In[25]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[26]:


save_data_to_db(data, paper_pmid)


# In[ ]:




