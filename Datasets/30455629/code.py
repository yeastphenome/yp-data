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


paper_pmid = 30455629
paper_name = 'fruhmann_cullin_2018' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[8]:


original_data = pd.read_excel('raw_data/Table_1_The Impact of ESCRT on Aβ1-42 Induced Membrane Lesions in a Yeast Model for Alzheimer’s Disease.XLS', 
                              sheet_name='Feuil1', skiprows=2)


# In[9]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[10]:


original_data.head()


# In[11]:


original_data['orf'] = original_data['systematic name'].astype(str)


# In[12]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[13]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[14]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[15]:


original_data['data'] = 0


# In[16]:


# Enhancers of AB-42 toxicity = more severe growth defect
ix = original_data['Growth'] <= 2
original_data.loc[ix,'data'] = original_data.loc[ix,'Growth'] - 3

# Suppressors of AB-42 toxicity = less severe growth defect
ix = original_data['Growth'] >= 3
original_data.loc[ix,'data'] = original_data.loc[ix,'Growth'] - 2


# In[17]:


original_data.set_index('orf', inplace=True)


# # Prepare the final dataset

# In[28]:


dataset_ids = [16614]


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


# In[38]:


save_data_to_db(data, paper_pmid)


# In[ ]:




