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


paper_pmid = 31647103
paper_name = 'schmidt_hombauer_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[6]:


original_data1 = pd.read_csv('raw_data/CAN1.txt', header=None, names=['genes'])
original_data2 = pd.read_csv('raw_data/lys2-10A.txt', header=None, names=['genes'])


# In[7]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[8]:


original_data1['genes'] = original_data1['genes'].astype(str)
original_data2['genes'] = original_data2['genes'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data1['genes'] = clean_genename(original_data1['genes'])
original_data2['genes'] = clean_genename(original_data2['genes'])


# In[10]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['genes'], to='orf')
original_data2['orfs'] = translate_sc(original_data2['genes'], to='orf')


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[13]:


original_data1['data'] = 1
original_data2['data'] = 1


# In[14]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)


# # Load & process tested strains

# In[15]:


tested = pd.read_excel('raw_data/transomic collection.xlsx', sheet_name='list with names')


# In[16]:


tested['SystematicName'] = tested['SystematicName'].astype(str)


# In[17]:


tested['SystematicName'] = clean_orf(tested['SystematicName'])


# In[18]:


tested['SystematicName'] = translate_sc(tested['SystematicName'], to='orf')


# In[19]:


t = looks_like_orf(tested['SystematicName'])
print(tested.loc[~t,])


# In[20]:


tested_orfs = tested['SystematicName'].unique()


# In[22]:


missing = [orf for orf in original_data1.index.values if orf not in tested_orfs]
missing


# In[23]:


missing = [orf for orf in original_data2.index.values if orf not in tested_orfs]
missing


# # Prepare the final dataset

# In[24]:


dataset_ids = [16441, 16442]


# In[25]:


datasets = datasets.reindex(index=dataset_ids)


# In[29]:


data = pd.DataFrame(index=tested_orfs, columns=datasets['name'].values, data=0)


# In[30]:


data.head()


# In[31]:


data.loc[original_data1.index, datasets['name'].values[0]] = original_data1['data']


# In[32]:


data.loc[original_data2.index, datasets['name'].values[1]] = original_data2['data']


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




