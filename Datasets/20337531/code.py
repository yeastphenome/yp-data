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


paper_pmid = 20337531
paper_name = 'dias_sa_correia_2010' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data1 = pd.read_csv('raw_data/TableS1.txt', header=None, names=['genes','orfs'], sep='\t')
original_data2 = pd.read_csv('raw_data/TableS2.txt', header=None, names=['genes','orfs'], sep='\t')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[7]:


original_data1['orfs'] = original_data1['orfs'].astype(str)
original_data2['orfs'] = original_data2['orfs'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data1['orfs'] = clean_orf(original_data1['orfs'])
original_data2['orfs'] = clean_orf(original_data2['orfs'])


# In[9]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['orfs'], to='orf')
original_data2['orfs'] = translate_sc(original_data2['orfs'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[12]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)


# In[13]:


# High-susceptibility
original_data1['data'] = -2

# Moderate-susceptibility
original_data2['data'] = -1


# In[14]:


original_data = original_data1.join(original_data2, lsuffix='_high', rsuffix='_mod', how='outer')


# In[15]:


original_data['data'] = original_data[['data_high','data_mod']].mean(axis=1)


# # Load & process tested strains

# In[16]:


tested = pd.read_excel('raw_data/List of strains tested.xlsx', sheet_name='Tabelle2')


# In[17]:


tested['ORF'] = clean_orf(tested['ORF'])


# In[18]:


tested = tested['ORF'].unique()


# In[19]:


tested = translate_sc(tested, to='orf')


# In[20]:


# Make sure everything translated ok
t = looks_like_orf(tested)
print(np.array(tested)[~np.array(t),])


# In[21]:


missing = [orf for orf in original_data.index.values if orf not in tested]


# In[22]:


# Remove missing (the data list contains HAP and HET screen results, so some of the genes are likely to be essential)
original_data.drop(index=missing, inplace=True)


# In[23]:


original_data.shape


# In[28]:


original_data = original_data.groupby(original_data.index).mean()


# In[29]:


original_data.shape


# # Prepare the final dataset

# In[24]:


dataset_ids = [16607]


# In[25]:


datasets = datasets.reindex(index=dataset_ids)


# In[26]:


data = pd.DataFrame(index=tested, columns=datasets['name'].values, data=0)


# In[30]:


data.loc[original_data.index, datasets['name'].values[0]] = original_data['data']


# In[31]:


data = data.groupby(data.index).mean()


# In[32]:


# Create row index
data.index.name='orf'


# In[33]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[40]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[41]:


from IO.save_data_to_db2 import *


# In[42]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[43]:


save_data_to_db(data, paper_pmid)


# In[ ]:




