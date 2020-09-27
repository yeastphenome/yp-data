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


paper_pmid = 15509654
paper_name = 'perrone_dawes_2005' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[18]:


original_data = pd.read_excel('raw_data/supp_tables.xlsx', sheet_name='Sheet1')


# In[19]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[20]:


original_data['orf'] = original_data['orf'].astype(str)


# In[21]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[22]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[23]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[24]:


original_data = original_data.loc[t,:]


# In[27]:


original_data = original_data.loc[original_data['fold'].notnull(),:]


# In[28]:


# Transforming into log-scale, so that WT (and all other mutants) have a value of 0 (fold = 1).
original_data['data'] = np.log2(original_data['fold'].values)


# In[29]:


original_data.set_index('orf', inplace=True)


# # Prepare the final dataset

# In[31]:


dataset_ids = [16434]


# In[32]:


datasets = datasets.reindex(index=dataset_ids)


# In[34]:


data = original_data[['data']].copy()


# In[35]:


data.columns = datasets['name'].values


# In[36]:


data = data.groupby(data.index).mean()


# In[37]:


# Create row index
data.index.name='orf'


# In[38]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[39]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[40]:


from IO.save_data_to_db2 import *


# In[41]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[42]:


save_data_to_db(data, paper_pmid)


# In[ ]:




