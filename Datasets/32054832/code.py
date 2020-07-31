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


paper_pmid = 32054832
paper_name = 'chen_petranovic_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data1 = pd.read_excel('raw_data/complete score data for SGA.xlsx', sheet_name='nonessential repeat 1')
original_data2 = pd.read_excel('raw_data/complete score data for SGA.xlsx', sheet_name='nonessential repeat 2')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[7]:


original_data1['orfs'] = original_data1['Array ORF'].astype(str)
original_data2['orfs'] = original_data2['Array ORF'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data1['orfs'] = clean_orf(original_data1['orfs'])
original_data2['orfs'] = clean_orf(original_data2['orfs'])


# In[9]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['orfs'], to='orf')
original_data2['orfs'] = translate_sc(original_data2['orfs'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[13]:


data1 = original_data1[['orfs','Score']].groupby('orfs').mean()
data2 = original_data2[['orfs','Score']].groupby('orfs').mean()


# In[14]:


data = data1.join(data2, lsuffix='_1', rsuffix='_2')


# In[15]:


data['data'] = data[['Score_1','Score_2']].mean(axis=1)


# In[18]:


data = data.drop(columns=['Score_1','Score_2'])


# # Prepare the final dataset

# In[19]:


dataset_ids = [16408]


# In[20]:


datasets = datasets.reindex(index=dataset_ids)


# In[21]:


data.columns = datasets['name'].values


# In[22]:


# Create row index
data.index.name='orf'


# In[23]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[24]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[25]:


from IO.save_data_to_db2 import *


# In[26]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[27]:


save_data_to_db(data, paper_pmid)


# In[ ]:




