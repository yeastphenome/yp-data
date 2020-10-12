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


paper_pmid = 20963216
paper_name = 'ratnakumar_oliver_2011' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/c0mb00114g.xls', sheet_name='A5', skiprows=1)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


orf_col = 'Systematic ID'


# In[8]:


original_data[orf_col] = original_data[orf_col].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data[orf_col] = clean_orf(original_data[orf_col])


# In[10]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data[orf_col], to='orf')


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[13]:


# Transform into TRT/UNT ratio, and take the log2
original_data['data'] = np.log2(-1/original_data['Fold-change'])


# In[14]:


original_data.set_index('orfs', inplace=True)


# # Prepare the final dataset

# In[19]:


dataset_ids = [16474]


# In[20]:


datasets = datasets.reindex(index=dataset_ids)


# In[21]:


data = original_data[['data']].copy()


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




