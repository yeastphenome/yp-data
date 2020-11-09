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


paper_pmid = 25074250
paper_name = 'hwang_naganuma_2014' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[6]:


original_data1 = pd.read_excel('raw_data/tables_1_2.xlsx', sheet_name='Resistance')
original_data2 = pd.read_excel('raw_data/tables_1_2.xlsx', sheet_name='Sensitivity')


# In[7]:


original_data = pd.concat([original_data1, original_data2], axis=0)


# In[8]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[11]:


original_data['orfs'] = original_data['ORF'].astype(str)


# In[12]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[15]:


original_data = original_data.reset_index()


# In[16]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[17]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[19]:


# Subtracting 1, so that the non-hit data can be automatically assigned to 0 at future analysis steps
original_data['data'] = original_data['IC50 of deletion cells/IC50 of control cells'] - 1


# In[20]:


original_data.set_index('orfs', inplace=True)


# # Prepare the final dataset

# In[21]:


dataset_ids = [16508]


# In[22]:


datasets = datasets.reindex(index=dataset_ids)


# In[24]:


data = original_data[['data']].copy()


# In[25]:


data.columns = datasets['name'].values


# In[26]:


data = data.groupby(data.index).mean()


# In[27]:


# Create row index
data.index.name='orf'


# In[28]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[29]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[30]:


from IO.save_data_to_db2 import *


# In[31]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[32]:


save_data_to_db(data, paper_pmid)


# In[ ]:




