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


paper_pmid = 24151994
paper_name = 'marek_korona_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# ### Maximum growth rate (MGR)

# In[7]:


original_data1 = pd.read_excel('raw_data/evo12196-sup-0001-tables1.xls', sheet_name='Table S5', skiprows=1)


# In[8]:


print('Original data dimensions: %d x %d' % (original_data1.shape))


# In[10]:


original_data1['ORF'] = original_data1['ORF'].astype(str)


# In[11]:


# Eliminate all white spaces & capitalize
original_data1['ORF'] = clean_orf(original_data1['ORF'])


# In[12]:


# Translate to ORFs 
original_data1['ORF'] = translate_sc(original_data1['ORF'], to='orf')


# In[14]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['ORF'])
print(original_data1.loc[~t,])


# In[16]:


original_data1['data'] = original_data1[['Replica 1.1','Replica 2.1','Replica 3.1']].mean(axis=1)


# In[17]:


original_data1.set_index('ORF', inplace=True)


# ### Maximum chronological lifespan (MLS)

# In[19]:


original_data2 = pd.read_excel('raw_data/evo12196-sup-0003-tables3.xls', sheet_name='3678 del', skiprows=1)


# In[20]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[21]:


original_data2['ORF'] = original_data2['ORF'].astype(str)


# In[22]:


# Eliminate all white spaces & capitalize
original_data2['ORF'] = clean_orf(original_data2['ORF'])


# In[23]:


# Translate to ORFs 
original_data2['ORF'] = translate_sc(original_data2['ORF'], to='orf')


# In[24]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['ORF'])
print(original_data2.loc[~t,])


# In[26]:


original_data2['data'] = original_data2[['Replica 1.1','Replica 2.1','Replica 3.1']].mean(axis=1)


# In[27]:


original_data2.set_index('ORF', inplace=True)


# ### Combine

# In[28]:


data = original_data1[['data']].join(original_data2[['data']], how='outer', lsuffix='_mgr', rsuffix='_mls')


# In[29]:


data.head()


# # Prepare the final dataset

# In[30]:


dataset_ids = [16564, 16565]


# In[31]:


datasets = datasets.reindex(index=dataset_ids)


# In[32]:


data.columns = datasets['name'].values


# In[33]:


data = data.groupby(data.index).mean()


# In[34]:


# Create row index
data.index.name='orf'


# In[35]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[37]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[38]:


from IO.save_data_to_db2 import *


# In[39]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[40]:


save_data_to_db(data, paper_pmid)


# In[ ]:




