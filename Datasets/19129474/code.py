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


# In[2]:


get_ipython().system('pwd')


# # Initial setup

# In[3]:


paper_pmid = 19129474
paper_name = 'tan_perrone_2009' 


# In[4]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[5]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[12]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, names=['genes'])


# In[13]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[14]:


original_data['genes'] = original_data['genes'].astype(str)


# In[15]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[16]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[17]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[18]:


original_data['data'] = -1


# In[19]:


original_data.set_index('orfs', inplace=True)


# # Load & process tested strains

# In[20]:


# Not available


# # Prepare the final dataset

# In[30]:


dataset_ids = [16450]


# In[31]:


datasets = datasets.reindex(index=dataset_ids)


# In[34]:


data = original_data['data'].to_frame().copy()
data.columns = dataset_ids


# In[35]:


data = data.groupby(data.index).mean()


# In[37]:


data.columns = datasets.loc[data.columns, 'name']


# In[38]:


# Create row index
data.index.name='orf'


# In[39]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[40]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[41]:


from IO.save_data_to_db2 import *


# In[43]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[44]:


save_data_to_db(data, paper_pmid)


# In[ ]:




