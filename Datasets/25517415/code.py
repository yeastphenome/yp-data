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


paper_pmid = 25517415
paper_name = 'winter_curtin_2014' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/Table_S1.xlsx', sheet_name='Sheet1')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data['orfs'] = original_data['ORF'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_genename(original_data['orfs'])


# In[9]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[11]:


data_dict = {'LOW': -2, 'MODERATELY LOW': -1, 'MODERATELY HIGH': 1, 'HIGH': 2}


# In[12]:


original_data['data'] = original_data['PHENOTYPE'].apply(lambda x: data_dict[x])


# In[13]:


original_data.set_index('orfs', inplace=True)


# # Prepare the final dataset

# In[15]:


dataset_ids = [16551]


# In[16]:


datasets = datasets.reindex(index=dataset_ids)


# In[17]:


data = original_data[['data']].copy()


# In[18]:


data.columns = datasets['name'].values


# In[19]:


data = data.groupby(data.index).mean()


# In[20]:


# Create row index
data.index.name='orf'


# In[21]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[22]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[23]:


from IO.save_data_to_db2 import *


# In[24]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[25]:


save_data_to_db(data, paper_pmid)


# In[ ]:




