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


# In[2]:


get_ipython().system('pwd')


# # Initial setup

# In[3]:


paper_pmid = 26994103
paper_name = 'luo_jiang_2016' 


# In[4]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[5]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[6]:


original_data1 = pd.read_excel('raw_data/Supplementary Table 2.xlsx', sheet_name='Sheet1', skiprows=2, names=['orf','gene','c1','c2'])
original_data2 = pd.read_excel('raw_data/Supplementary Table 3.xlsx', sheet_name='Sheet1', skiprows=2, names=['orf','gene','c1'])


# In[7]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[8]:


original_data1['orf'] = original_data1['orf'].astype(str)
original_data2['orf'] = original_data2['orf'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data1['orf'] = clean_orf(original_data1['orf'])
original_data2['orf'] = clean_orf(original_data2['orf'])


# In[10]:


# Translate to ORFs 
original_data1['orf'] = translate_sc(original_data1['orf'], to='orf')
original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')


# In[11]:


# Make sure everything translated ok
t1 = looks_like_orf(original_data1['orf'])
t2 = looks_like_orf(original_data2['orf'])


# In[12]:


print(original_data1.loc[~t1,])


# In[13]:


print(original_data2.loc[~t2,])


# In[14]:


original_data1.set_index('orf', inplace=True)
original_data2.set_index('orf', inplace=True)


# In[15]:


original_data1['data'] = -1
original_data2['data'] = 1


# In[16]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[17]:


original_data = original_data[['data_1','data_2']].copy()


# In[18]:


original_data[original_data.isnull()] = 0


# In[20]:


dataset_ids = [16448,16449]
datasets = datasets.reindex(index=dataset_ids)


# In[21]:


original_data.columns = dataset_ids


# In[22]:


print('Final data dimensions: %d x %d' % (original_data.shape))


# # Load the tested strains

# In[19]:


tested = pd.read_excel('raw_data/DELETION LIBRARY.xlsx', sheet_name='DELETION LIBRARY', skiprows=1)
tested = np.array(tested['ORF name'].unique())
tested = clean_orf(tested)
tested = translate_sc(tested, to='orf')
tested[tested == 'YELOO1C'] = 'YEL001C'
# Make sure everything translated ok
t = looks_like_orf(tested)
print(tested[~np.array(t)])
tested = np.setdiff1d(tested, np.array(['YMR41W']))


# # Prepare the final dataset

# In[20]:





# In[21]:


data = pd.DataFrame(index=tested, columns=datasets['name'].values, data=0)


# In[22]:


data.loc[original_data.index, datasets['name'].values[0]] = original_data['data_1']


# In[23]:


data.loc[original_data.index, datasets['name'].values[1]] = original_data['data_2']


# In[24]:


# Create row index
data.index.name='orf'


# In[25]:


# If the same strain is present more than once, average its values
data = data.groupby('orf').mean()


# In[26]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[27]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[28]:


from IO.save_data_to_db2 import *


# In[29]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[30]:


save_data_to_db(data, paper_pmid)

