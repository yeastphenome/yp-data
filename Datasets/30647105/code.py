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


paper_pmid = 30647105
paper_name = 'alhoch_tang_2019' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[17]:


original_data = pd.read_excel('raw_data/S6 SC genomic screen BHA and BPA TN1.xlsx', sheet_name='Sc BPA and BHA genomic', skiprows=1)


# In[18]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[19]:


original_data['ORF name'] = original_data['ORF name'].astype(str)


# In[20]:


# Eliminate all white spaces & capitalize
original_data['ORF name'] = clean_orf(original_data['ORF name'])


# In[21]:


# Translate to ORFs 
original_data['ORF name'] = translate_sc(original_data['ORF name'], to='orf')


# In[22]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF name'])
print(original_data.loc[~t,])


# In[23]:


original_data = original_data.loc[t,]


# In[24]:


original_data.set_index('ORF name', inplace=True)


# In[25]:


original_data[original_data.notnull()] = -1


# In[26]:


original_data[original_data.isnull()] = 0


# In[27]:


original_data.sum(axis=0)


# In[33]:


original_data = original_data.astype(float)


# # Prepare the final dataset

# In[34]:


dataset_ids = [16600, 16599]


# In[35]:


datasets = datasets.reindex(index=dataset_ids)


# In[36]:


data = original_data.copy()


# In[37]:


data.columns = datasets['name'].values


# In[38]:


data = data.groupby(data.index).mean()


# In[39]:


# Create row index
data.index.name='orf'


# In[40]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[42]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[43]:


from IO.save_data_to_db2 import *


# In[44]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[45]:


save_data_to_db(data, paper_pmid)


# In[ ]:




