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


paper_pmid = 18975972
paper_name = 'mccue_phang_2008' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[6]:


original_data = pd.read_csv('raw_data/hits.txt', header=0, sep='\t')


# In[7]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[9]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[10]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[11]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'], to='orf')


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[13]:


original_data['data'] = -1


# In[14]:


original_data.set_index('ORF', inplace=True)


# # Prepare the final dataset

# In[16]:


dataset_ids = [16437]


# In[17]:


datasets = datasets.reindex(index=dataset_ids)


# In[19]:


data = original_data[['data']].copy()


# In[20]:


data.columns = datasets['name'].values


# In[21]:


data = data.groupby(data.index).mean()


# In[22]:


# Create row index
data.index.name='orf'


# In[23]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[24]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[25]:


from IO.save_data_to_db2 import *


# In[26]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[27]:


save_data_to_db(data, paper_pmid)


# In[ ]:




