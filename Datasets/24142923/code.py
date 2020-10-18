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


paper_pmid = 24142923
paper_name = 'jarolim_dawes_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[12]:


original_data1 = pd.read_excel('raw_data/FileS1.xlsx', sheet_name='Exponential_Resistant', skiprows=2)
original_data2 = pd.read_excel('raw_data/FileS1.xlsx', sheet_name='Exponential_Sensitive', skiprows=2)
original_data3 = pd.read_excel('raw_data/FileS1.xlsx', sheet_name='Stationary_Resistant', skiprows=2)
original_data4 = pd.read_excel('raw_data/FileS1.xlsx', sheet_name='Stationary_Sensitive', skiprows=2)


# In[13]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))
print('Original data dimensions: %d x %d' % (original_data3.shape))
print('Original data dimensions: %d x %d' % (original_data4.shape))


# In[14]:


orf_col = 'Systematic Name'
original_data1['orfs'] = original_data1['SystematicName'].astype(str)
original_data2['orfs'] = original_data2[orf_col].astype(str)
original_data3['orfs'] = original_data3[orf_col].astype(str)
original_data4['orfs'] = original_data4[orf_col].astype(str)


# In[15]:


# Eliminate all white spaces & capitalize
original_data1['orfs'] = clean_orf(original_data1['orfs'])
original_data2['orfs'] = clean_orf(original_data2['orfs'])
original_data3['orfs'] = clean_orf(original_data3['orfs'])
original_data4['orfs'] = clean_orf(original_data4['orfs'])


# In[16]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['orfs'], to='orf')
original_data2['orfs'] = translate_sc(original_data2['orfs'], to='orf')
original_data3['orfs'] = translate_sc(original_data3['orfs'], to='orf')
original_data4['orfs'] = translate_sc(original_data4['orfs'], to='orf')


# In[20]:


# Fix typos
original_data3.loc[original_data3['orfs']=='YML095-A','orfs'] = 'YML095C-A'


# In[17]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[18]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[21]:


# Make sure everything translated ok
t = looks_like_orf(original_data3['orfs'])
print(original_data3.loc[~t,])


# In[22]:


# Make sure everything translated ok
t = looks_like_orf(original_data4['orfs'])
print(original_data4.loc[~t,])


# In[23]:


original_data1.head()


# In[24]:


original_data1['data'] = original_data1['Rating']
original_data2['data'] = original_data2['Rating']
original_data3['data'] = original_data3['Rating']
original_data4['data'] = original_data4['Rating']


# In[25]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)
original_data3.set_index('orfs', inplace=True)
original_data4.set_index('orfs', inplace=True)


# In[26]:


data1 = original_data1[['data']].join(original_data2[['data']], how='outer', lsuffix='_r', rsuffix='_s')


# In[29]:


data1['data'] = data1[['data_s','data_r']].mean(axis=1)


# In[31]:


data2 = original_data3[['data']].join(original_data4[['data']], how='outer', lsuffix='_r', rsuffix='_s')


# In[32]:


data2['data'] = data2[['data_s','data_r']].mean(axis=1)


# In[34]:


data = data1[['data']].join(data2[['data']], how='outer', lsuffix='_exp', rsuffix='_stn')


# In[37]:


data[data.isnull()] = 0


# In[40]:


data.head()


# # Prepare the final dataset

# In[41]:


dataset_ids = [16624,16538]


# In[42]:


datasets = datasets.reindex(index=dataset_ids)


# In[43]:


data.columns = datasets['name'].values


# In[44]:


data = data.groupby(data.index).mean()


# In[45]:


# Create row index
data.index.name='orf'


# In[46]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[48]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[49]:


from IO.save_data_to_db2 import *


# In[50]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[51]:


save_data_to_db(data, paper_pmid)


# In[ ]:




