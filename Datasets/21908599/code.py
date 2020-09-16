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

# In[3]:


paper_pmid = 21908599
paper_name = 'minear_cyert_2011' 


# In[4]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[5]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[6]:


original_data = pd.read_excel('raw_data/TableS1.xlsx', sheet_name='Sheet1')


# In[7]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[8]:


original_data.head()


# In[9]:


original_data['Locus'] = original_data['Locus'].astype(str)


# In[10]:


# Eliminate all white spaces & capitalize
original_data['Locus'] = clean_orf(original_data['Locus'])


# In[11]:


# Translate to ORFs 
original_data['Locus'] = translate_sc(original_data['Locus'], to='orf')


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data['Locus'])
print(original_data.loc[~t,])


# In[13]:


original_data['data'] = -original_data['Defect vs WT']


# In[15]:


original_data.set_index('Locus', inplace=True)


# # Prepare the final dataset

# In[16]:


dataset_ids = [16484]


# In[17]:


datasets = datasets.reindex(index=dataset_ids)


# In[18]:


data = original_data['data'].to_frame()


# In[19]:


data.columns = datasets['name'].values


# In[20]:


data = data.groupby(data.index).mean()


# In[21]:


# Create row index
data.index.name='orf'


# In[22]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[23]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[24]:


from IO.save_data_to_db2 import *


# In[25]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[26]:


save_data_to_db(data, paper_pmid)


# In[ ]:




