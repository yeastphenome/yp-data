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


paper_pmid = 23103169
paper_name = 'qian_zhang_2012' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[7]:


original_data = pd.read_excel('raw_data/mmc2.xls', sheet_name='fitness_combined', skiprows=1)


# In[8]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[10]:


original_data = original_data[['ORF','YPD rep1 fitness','YPD rep2 fitness','YPG fitness','YPE fitness','SC fitness','OAK fitness','ETH fitness']]


# In[11]:


original_data.head()


# In[12]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[13]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[14]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'], to='orf')


# In[15]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[16]:


original_data['YPD fitness'] = original_data[['YPD rep1 fitness','YPD rep2 fitness']].mean(axis=1)


# In[20]:


original_data.drop(columns=['YPD rep1 fitness','YPD rep2 fitness'], inplace=True)


# In[18]:


original_data = original_data.groupby('ORF').mean()


# # Prepare the final dataset

# In[22]:


dataset_ids = [16489, 16488, 16490, 16492, 16491, 16487]


# In[23]:


datasets = datasets.reindex(index=dataset_ids)


# In[24]:


data = original_data.copy()


# In[25]:


data.columns = datasets['name'].values


# In[26]:


# Create row index
data.index.name='orf'


# In[27]:


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




