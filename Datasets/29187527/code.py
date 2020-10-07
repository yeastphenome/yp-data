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


paper_pmid = 29187527
paper_name = 'eisenberg_bord_bohnert_2018' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# ### Part 1

# In[24]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, names=['genes','data'], sep='\t')


# In[25]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[26]:


original_data['genes'] = original_data['genes'].astype(str)


# In[27]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[28]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[29]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[30]:


original_data.set_index('orfs', inplace=True)


# ### Part 2

# In[35]:


original_data2 = pd.read_excel('raw_data/JCB_201704122_TableS1.xlsx', sheet_name='HITS from screens', skiprows=2, nrows=30)


# In[38]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[39]:


original_data2['ORF'] = original_data2['ORF'].astype(str)


# In[40]:


# Eliminate all white spaces & capitalize
original_data2['ORF'] = clean_orf(original_data2['ORF'])


# In[41]:


# Translate to ORFs 
original_data2['ORF'] = translate_sc(original_data2['ORF'], to='orf')


# In[42]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['ORF'])
print(original_data2.loc[~t,])


# In[43]:


original_data2['data'] = 1


# In[44]:


original_data2.set_index('ORF', inplace=True)


# ### Merge

# In[54]:


data = original_data[['data']].join(original_data2[['data']], lsuffix='_lower', rsuffix='_higher', how='outer')


# In[55]:


data['data'] = data[['data_lower','data_higher']].mean(axis=1)


# In[57]:


data = data[['data']]


# # Prepare the final dataset

# In[58]:


dataset_ids = [16438]


# In[59]:


datasets = datasets.reindex(index=dataset_ids)


# In[60]:


data.columns = datasets['name'].values


# In[61]:


data = data.groupby(data.index).mean()


# In[62]:


# Create row index
data.index.name='orf'


# In[63]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[65]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[66]:


from IO.save_data_to_db2 import *


# In[67]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[68]:


save_data_to_db(data, paper_pmid)


# In[ ]:




