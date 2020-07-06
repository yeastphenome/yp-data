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


paper_pmid = 19416971
paper_name = 'khozoie_avery_2009' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[8]:


original_data = pd.read_excel('raw_data/jbc.M109.005843-1.xls', sheet_name='Initial QN Screen', skiprows=3)


# In[9]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[10]:


original_data.head()


# In[12]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[13]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[14]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'], to='orf')


# In[16]:


to_drop = original_data.loc[original_data['ORF']=='EMPTY',]


# In[18]:


original_data.drop(index=to_drop.index, inplace=True)


# In[19]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[20]:


original_data.drop(index=original_data.loc[~t,:].index, inplace=True)


# In[24]:


# Eliminate the slow growin strains for which no accurate growth ratio could be calculated
to_drop = original_data.loc[original_data['Adjusted GR']=='SLOW',]


# In[26]:


original_data.drop(index=to_drop.index, inplace=True)


# In[27]:


original_data.shape


# In[36]:


data = original_data[['ORF','Adjusted GR']].copy()


# In[37]:


data.set_index('ORF', inplace=True)


# In[38]:


data['Adjusted GR'] = data['Adjusted GR'].astype(float)


# # Prepare the final dataset

# In[39]:


dataset_ids = [16533]


# In[40]:


datasets = datasets.reindex(index=dataset_ids)


# In[41]:


data.columns = datasets['name'].values


# In[42]:


data = data.groupby(data.index).mean()


# In[43]:


# Create row index
data.index.name='orf'


# In[45]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[46]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[47]:


from IO.save_data_to_db2 import *


# In[48]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[49]:


save_data_to_db(data, paper_pmid)


# In[ ]:




