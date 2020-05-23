#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd

import sys

from os.path import expanduser
sys.path.append(expanduser('~') + '/Lab/Utils/Python/')

from Conversions.translate import *
from Strings.is_a import *


# # Initial setup

# In[3]:


paper_pmid = 31029968
paper_name = 'alfatah_arumugam_2019' 


# In[4]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[5]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[6]:


original_data = pd.read_excel('raw_data/Supplementary table. 1.xlsx', sheet_name='HOP_Monoethylhexylphthalicacid_')


# In[7]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[8]:


original_data['genes'] = original_data['genes'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[10]:


# If possible, fix typos, omissions, etc.
original_data.loc[original_data['genes'].str.contains('2001-10-01'),'genes'] = 'OCT1'
original_data.loc[original_data['genes'].str.contains('YBR160WAS'),'genes'] = 'YBR160W'


# In[11]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])


# In[13]:


print(original_data.loc[~t,])


# In[14]:


# If the same strain is present more than once, average its values
data = original_data.groupby('orfs')['logFC'].mean().to_frame()


# In[15]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Prepare the final dataset

# In[16]:


dataset_ids = [16439]


# In[17]:


datasets = datasets.reindex(index=dataset_ids)


# In[18]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[19]:


# Create row index
data.index.name='orf'


# # Print out

# In[20]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[21]:


from IO.save_data_to_db2 import *


# In[22]:


save_data_to_db(data, paper_pmid)

