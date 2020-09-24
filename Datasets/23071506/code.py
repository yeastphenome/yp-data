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

# In[3]:


paper_pmid = 23071506
paper_name = 'lockshon_kennedy_2012' 


# In[4]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[5]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[6]:


original_data = pd.read_excel('raw_data/journal.pone.0045049.s009.XLS', sheet_name='PO inhib order', skiprows=1)


# In[7]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[8]:


original_data['orf'] = original_data['orf'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[10]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[13]:


original_data['PO inh'] = -original_data['PO inh']


# In[16]:


((original_data['PO inh'] < 0) & (original_data['POstim'] > 0)).sum()


# In[14]:


original_data['data'] = original_data[['PO inh','POstim']].sum(axis=1)


# In[17]:


original_data.head()


# In[18]:


original_data.set_index('orf', inplace=True)


# # Prepare the final dataset

# In[19]:


dataset_ids = [16611]


# In[20]:


datasets = datasets.reindex(index=dataset_ids)


# In[21]:


data = original_data['data'].to_frame()


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




