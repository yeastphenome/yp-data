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


paper_pmid = 20960220
paper_name = 'jayakody_kitagaki_2011' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/10529_2010_437_MOESM1_ESM.xlsx', sheet_name='Table 1', skiprows=7)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[9]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[14]:


# Convert data to follow convention (centered on 0, negative = defect)
original_data['data'] = original_data['Ratio']-1


# In[12]:


original_data.set_index('ORF', inplace=True)


# # Prepare the final dataset

# In[15]:


dataset_ids = [16473]


# In[16]:


datasets = datasets.reindex(index=dataset_ids)


# In[17]:


data = original_data['data'].to_frame()


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




