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


paper_pmid = 19635755
paper_name = 'jo_vulpe_2009' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[10]:


original_data1 = pd.read_excel('raw_data/Table_S2.xlsx', sheet_name='Table 1', skiprows=4, header=None)
original_data2 = pd.read_excel('raw_data/Table_S3.xlsx', sheet_name='Table 1', skiprows=4, header=None)


# In[11]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[15]:


original_data1.head()


# In[16]:


original_data1.columns = ['orfs','genes','37uM_5gen','75uM_5gen','150uM_5gen','37uM_15gen','75uM_15gen','150uM_15gen', 'num_hits','t']
original_data2.columns = ['orfs','genes','75uM_5gen','150uM_5gen','300uM_5gen','75ug_15gen','150uM_15gen','300uM_15gen', 'num_hits','t']


# In[19]:


original_data1['orfs'] = original_data1['orfs'].astype(str)
original_data2['orfs'] = original_data2['orfs'].astype(str)


# In[20]:


# Eliminate all white spaces & capitalize
original_data1['orfs'] = clean_orf(original_data1['orfs'])
original_data2['orfs'] = clean_orf(original_data2['orfs'])


# In[21]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['orfs'], to='orf')
original_data2['orfs'] = translate_sc(original_data2['orfs'], to='orf')


# In[22]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[23]:


original_data1 = original_data1.loc[t,]


# In[24]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[25]:


original_data2 = original_data2.loc[t,]


# In[28]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)


# In[30]:


original_data = original_data1[['37uM_5gen','75uM_5gen','150uM_5gen','37uM_15gen','75uM_15gen','150uM_15gen']].join(
    original_data2[['75uM_5gen','150uM_5gen','300uM_5gen','75ug_15gen','150uM_15gen','300uM_15gen']], 
    lsuffix='_mma', rsuffix='_as', how='outer')


# In[33]:


original_data[original_data.isnull()] = 0


# In[41]:


original_data = original_data.astype(float)


# # Prepare the final dataset

# In[42]:


dataset_ids = [16656,16655,16654,16653,16652,16651,16659,16658,16657,16662,16661,16660]


# In[43]:


datasets = datasets.reindex(index=dataset_ids)


# In[44]:


data = original_data.copy()


# In[45]:


data.columns = datasets['name'].values


# In[46]:


data = data.groupby(data.index).mean()


# In[47]:


# Create row index
data.index.name='orf'


# In[48]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[49]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[50]:


from IO.save_data_to_db2 import *


# In[51]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[52]:


save_data_to_db(data, paper_pmid)


# In[ ]:




