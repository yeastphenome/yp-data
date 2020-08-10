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


paper_pmid = 24211263
paper_name = 'teng_hardwick_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[19]:


original_data = pd.read_excel('raw_data/mmc5.xlsx', sheet_name='749 YKOs with low aa overgrowth')


# In[20]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[21]:


original_data = original_data.rename(columns={'ORF name of YKOs with overgrowth phenotype':'orfs'})


# In[22]:


original_data['orfs'] = original_data['orfs'].astype(str)


# In[23]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[24]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[25]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[26]:


original_data = original_data.loc[t,:]


# In[27]:


original_data['data'] = 1


# In[28]:


original_data.set_index('orfs', inplace=True)


# # Prepare the final dataset

# In[39]:


dataset_ids = [16393]


# In[40]:


datasets = datasets.reindex(index=dataset_ids)


# In[51]:


data = original_data['data'].to_frame()


# In[52]:


data.columns = datasets['name'].values


# In[53]:


data = data.groupby(data.index).mean()


# In[54]:


# Create row index
data.index.name='orf'


# In[55]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[56]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[57]:


from IO.save_data_to_db2 import *


# In[58]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[59]:


save_data_to_db(data, paper_pmid)


# In[ ]:




