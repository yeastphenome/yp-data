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


paper_pmid = 31181115
paper_name = 'wong_arumugam_2019' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/journal.pone.0218189.s008.xlsx', sheet_name='Plan 1', skiprows=2)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data['genes'] = original_data['genes'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[11]:


original_data.loc[original_data['genes']=='2001-10-01000000','genes'] = 'OCT1'


# In[12]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[13]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[14]:


original_data = original_data.loc[t,:]


# In[15]:


original_data.set_index('orfs', inplace=True)


# # Prepare the final dataset

# In[16]:


dataset_ids = [16598]


# In[17]:


datasets = datasets.reindex(index=dataset_ids)


# In[18]:


data = original_data['logFC'].to_frame()


# In[19]:


data.columns = datasets['name'].values


# In[20]:


data = data.groupby(data.index).mean()


# In[21]:


# Create row index
data.index.name='orf'


# In[22]:


print('Final data dimensions: %d x %d' % (data.shape))


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




