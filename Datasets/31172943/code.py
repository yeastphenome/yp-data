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


paper_pmid = 31172943
paper_name = 'dederer_lamberg_2019' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data1 = pd.read_excel('raw_data/elife-45506-supp1-v2.xlsx', sheet_name='tFT-Pex15delta30')
original_data2 = pd.read_excel('raw_data/elife-45506-supp1-v2.xlsx', sheet_name='tFT-Tom5')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[7]:


original_data1['ORF'] = original_data1['ORF'].astype(str)
original_data2['ORF'] = original_data2['ORF'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data1['ORF'] = clean_orf(original_data1['ORF'])
original_data2['ORF'] = clean_orf(original_data2['ORF'])


# In[9]:


# Translate to ORFs 
original_data1['ORF'] = translate_sc(original_data1['ORF'], to='orf')
original_data2['ORF'] = translate_sc(original_data2['ORF'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['ORF'])
print(original_data1.loc[~t,])


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['ORF'])
print(original_data2.loc[~t,])


# In[12]:


original_data1.set_index('ORF', inplace=True)
original_data2.set_index('ORF', inplace=True)


# In[20]:


original_data1 = original_data1.groupby(original_data1.index).mean()


# In[21]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[26]:


data = original_data1[['ratio.mean']].join(original_data2[['ratio.mean']], how='outer', lsuffix='_1', rsuffix='_2')


# # Prepare the final dataset

# In[29]:


dataset_ids = [16549,16550]


# In[30]:


datasets = datasets.reindex(index=dataset_ids)


# In[31]:


data.columns = datasets['name'].values


# In[33]:


data = data.groupby(data.index).mean()


# In[34]:


# Create row index
data.index.name='orf'


# In[35]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[36]:


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




