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


paper_pmid = 25359478
paper_name = 'jayakody_kitagaki_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/Final results 1 fileLahiruscreening.xls', sheet_name='All screen data', skiprows=8)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[9]:


original_data.columns


# In[10]:


original_data['orfs'] = original_data['Mutant '].astype(str)


# In[11]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[12]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[13]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[14]:


original_data['data'] = original_data['SI AVG']


# In[15]:


original_data.set_index('orfs', inplace=True)


# In[26]:


original_data['data'] = pd.to_numeric(original_data['data'], errors='coerce', downcast='float')


# # Prepare the final dataset

# In[28]:


dataset_ids = [16610]


# In[29]:


datasets = datasets.reindex(index=dataset_ids)


# In[30]:


data = original_data[['data']].copy()


# In[31]:


data.columns = datasets['name'].values


# In[32]:


data = data.groupby(data.index).mean()


# In[33]:


# Create row index
data.index.name='orf'


# In[34]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[35]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[36]:


from IO.save_data_to_db2 import *


# In[37]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[40]:


save_data_to_db(data, paper_pmid)


# In[ ]:




