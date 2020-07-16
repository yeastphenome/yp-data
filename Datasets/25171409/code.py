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


paper_pmid = 25171409
paper_name = 'jarosz_lindquist_2014' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data1 = pd.read_excel('raw_data/mmc3.xls', sheet_name='Sheet1')
original_data2 = pd.read_excel('raw_data/mmc4.xlsx', sheet_name='Sheet1')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[7]:


original_data1['ORF name'] = original_data1['ORF name'].astype(str)
original_data2['ORF name'] = original_data2['ORF name'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data1['ORF name'] = clean_orf(original_data1['ORF name'])
original_data2['ORF name'] = clean_orf(original_data2['ORF name'])


# In[9]:


# Translate to ORFs 
original_data1['ORF name'] = translate_sc(original_data1['ORF name'], to='orf')
original_data2['ORF name'] = translate_sc(original_data2['ORF name'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['ORF name'])
print(original_data1.loc[~t,])


# In[12]:


original_data1.drop(index=original_data1.loc[~t].index, inplace=True)


# In[13]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['ORF name'])
print(original_data2.loc[~t,])


# In[14]:


original_data2.drop(index=original_data2.loc[~t].index, inplace=True)


# In[15]:


original_data1['data'] = -1
original_data2['data'] = 1


# In[16]:


original_data1.set_index('ORF name', inplace=True)
original_data2.set_index('ORF name', inplace=True)


# In[21]:


data = original_data1['data'].to_frame().join(original_data2['data'].to_frame(), how='outer', lsuffix='_reduced', rsuffix='_enhanced')


# In[23]:


data.shape


# In[24]:


data['data'] = data[['data_reduced','data_enhanced']].mean(axis=1)


# In[31]:


data.drop(columns=['data_enhanced','data_reduced'], inplace=True)


# # Prepare the final dataset

# In[28]:


dataset_ids = [16510]


# In[29]:


datasets = datasets.reindex(index=dataset_ids)


# In[32]:


data.columns = datasets['name'].values


# In[33]:


data = data.groupby(data.index).mean()


# In[34]:


# Create row index
data.index.name='orf'


# In[35]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[36]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[37]:


from IO.save_data_to_db2 import *


# In[38]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[39]:


save_data_to_db(data, paper_pmid)


# In[ ]:




