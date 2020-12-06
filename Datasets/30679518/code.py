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


paper_pmid = 30679518
paper_name = 'alfatah_arumugam_2019' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/41598_2018_35979_MOESM2_ESM.xlsx', sheet_name='HOP_Hypoculoside')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data['genes'] = original_data['genes'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[9]:


# # If possible, fix typos, omissions, etc.
# original_data.loc[original_data['genes'].str.contains('2001-10-01'),'genes'] = 'OCT1'
# original_data.loc[original_data['genes'].str.contains('YBR160WAS'),'genes'] = 'YBR160W'


# In[10]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])


# In[12]:


print(original_data.loc[~t,])


# In[13]:


to_remove = original_data.loc[original_data['genes'] == '37165',:].index.values
original_data.drop(index=to_remove, inplace=True)


# In[14]:


original_data.shape


# In[22]:


# If the same strain is present more than once, average its values
data = original_data.groupby('orfs')['logFC'].mean().to_frame()


# In[23]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Prepare the final dataset

# In[24]:


dataset_ids = [16440]


# In[25]:


datasets = datasets.reindex(index=dataset_ids)


# In[27]:


data.columns = datasets['name'].values


# In[28]:


data = data.groupby(data.index).mean()


# In[29]:


# Create row index
data.index.name='orf'


# In[30]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[31]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[21]:


from IO.save_data_to_db2 import *


# In[ ]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[22]:


save_data_to_db(data, paper_pmid)

