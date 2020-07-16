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

# In[5]:


original_data = pd.read_excel('raw_data/jbc.M109.005843-1.xls', sheet_name='Initial QN Screen', skiprows=3)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[10]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'], to='orf')


# In[11]:


to_drop = original_data.loc[original_data['ORF']=='EMPTY',]


# In[12]:


original_data.drop(index=to_drop.index, inplace=True)


# In[13]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[14]:


original_data.drop(index=original_data.loc[~t,:].index, inplace=True)


# In[15]:


# Eliminate the slow growin strains for which no accurate growth ratio could be calculated
to_drop = original_data.loc[original_data['Adjusted GR']=='SLOW',]


# In[16]:


original_data.drop(index=to_drop.index, inplace=True)


# In[17]:


original_data.shape


# In[25]:


# Reverse the growth ratio so that lower values correspond to decreased growth and viceversa (originally, GR is reported as untreated vs treated)
original_data['GR2'] = 1 / original_data['Growth Ratio (GR)']


# In[26]:


# Normalize by plate median (as done oridinally)
def normalize_by_plate_median(plate_data):
    plate_median = plate_data['GR2'].median()
    plate_data['GR2_adjusted'] = plate_data['GR2'] / plate_median
    return plate_data


# In[27]:


original_data2 = original_data.groupby('Plate').apply(normalize_by_plate_median)


# In[28]:


original_data2.head()


# In[40]:


data = original_data2[['ORF','GR2_adjusted']].copy()


# In[41]:


data['GR2_adjusted'] = data['GR2_adjusted'].astype(float)


# In[42]:


data.set_index('ORF', inplace=True)


# # Prepare the final dataset

# In[43]:


dataset_ids = [16533]


# In[44]:


datasets = datasets.reindex(index=dataset_ids)


# In[45]:


data.columns = datasets['name'].values


# In[47]:


data = data.groupby(data.index).mean()


# In[48]:


# Create row index
data.index.name='orf'


# In[49]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[50]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[51]:


from IO.save_data_to_db2 import *


# In[52]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[53]:


save_data_to_db(data, paper_pmid)


# In[ ]:




