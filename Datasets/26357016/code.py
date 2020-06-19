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

# In[30]:


paper_pmid = 26357016
paper_name = 'frohlich_walther_2015' 


# In[31]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[32]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[36]:


original_data = pd.read_excel('raw_data/Myriocin screen.xlsx', sheet_name='Myriocin screen')


# In[37]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[38]:


# Keep only deletions
original_data = original_data.loc[original_data['Mutation']=='DELETION',:]


# In[39]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[40]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[41]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'].values, to='orf')


# In[42]:


original_data.loc[original_data['ORF']=='YLR287-A','ORF'] = 'YLR287C-A'


# In[43]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[47]:


data = original_data.groupby('ORF').mean()


# In[51]:


data = data.loc[:, ['WT.mean','MYR.mean']]


# In[52]:


# Normalize by untreated
data['MYR.mean'] = data['MYR.mean']/data['WT.mean']


# In[53]:


data.head()


# # Prepare the final dataset

# In[54]:


dataset_ids = [16496,16458]


# In[55]:


datasets = datasets.reindex(index=dataset_ids)


# In[56]:


data.columns = datasets['name'].values


# In[58]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[59]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[60]:


from IO.save_data_to_db2 import *


# In[61]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[62]:


save_data_to_db(data, paper_pmid)


# In[ ]:




