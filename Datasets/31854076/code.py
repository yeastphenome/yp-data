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


paper_pmid = 31854076
paper_name = 'novarina_chang_2019' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[29]:


original_data = pd.read_excel('raw_data/acel13084-sup-0001-files1.xlsx', sheet_nam='Filtered', header=8, nrows=4898-9)


# In[30]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[31]:


original_data['Systematic_Name'] = original_data['Systematic_Name'].astype(str)


# In[32]:


# Eliminate all white spaces & capitalize
original_data['Systematic_Name'] = clean_orf(original_data['Systematic_Name'])


# In[33]:


original_data.loc[original_data['Systematic_Name']=='YLR287-A','Systematic_Name'] = 'YLR287C-A'


# In[34]:


# Make sure everything translated ok
t = looks_like_orf(original_data['Systematic_Name'])


# In[26]:


print(original_data.loc[~t,])


# In[36]:


original_data.drop(index=original_data.loc[~t,].index, inplace=True)


# In[37]:


print(original_data.shape)


# In[38]:


original_data.head()


# In[39]:


original_data.set_index('Systematic_Name', inplace=True)


# In[40]:


# Get the relevant columns
original_data = original_data.loc[:, 'percentage_escapers'].to_frame()


# In[41]:


print(original_data.shape)


# In[42]:


# If the same strain is present more than once, average its values
data = original_data.groupby(original_data.index.values).mean()


# In[43]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Prepare the final dataset

# In[44]:


dataset_ids = [16401]


# In[45]:


datasets = datasets.reindex(index=dataset_ids)


# In[46]:


# Create row index
data.index.name='orf'
data.columns = datasets['name']


# # Print out

# In[48]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[49]:


from IO.save_data_to_db2 import *


# In[50]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[51]:


save_data_to_db(data, paper_pmid)


# In[ ]:




