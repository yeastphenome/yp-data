#!/usr/bin/env python
# coding: utf-8

# In[3]:


import numpy as np
import pandas as pd

import sys

from os.path import expanduser
sys.path.append(expanduser('~') + '/Lab/Utils/Python/')

from Conversions.translate import *
from Strings.is_a import *


# In[4]:


get_ipython().system('pwd')


# # Initial setup

# In[5]:


paper_pmid = 23157175
paper_name = 'zhang_jiang_2013' 


# In[6]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[7]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[8]:


original_data1 = pd.read_excel('raw_data/hits.xlsx', sheet_name='DMSO4')
original_data2 = pd.read_excel('raw_data/hits.xlsx', sheet_name='DMSO8')


# In[9]:


original_data1.columns = [x.strip() for x in original_data1.columns]
original_data2.columns = [x.strip() for x in original_data2.columns]


# In[11]:


original_data1.drop(index=original_data1.loc[original_data1['DMSO sensitivity'].isnull()].index, inplace=True)
original_data2.drop(index=original_data2.loc[original_data2['8% DMSO sensitivity'].isnull()].index, inplace=True)


# In[12]:


original_data1['DMSO sensitivity'] = original_data1['DMSO sensitivity'].apply(lambda x: len(x.strip()))
original_data2['8% DMSO sensitivity'] = original_data2['8% DMSO sensitivity'].apply(lambda x: len(x.strip()))


# In[13]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[14]:


original_data1['genes'] = original_data1['Gene'].astype(str)
original_data2['genes'] = original_data2['Gene'].astype(str)


# In[15]:


# Eliminate all white spaces & capitalize
original_data1['genes'] = clean_genename(original_data1['genes'])
original_data2['genes'] = clean_genename(original_data2['genes'])


# In[16]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['genes'].values, to='orf')
original_data2['orfs'] = translate_sc(original_data2['genes'].values, to='orf')


# In[17]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[18]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[19]:


original_data1['data'] = -original_data1['DMSO sensitivity']
original_data2['data'] = -original_data2['8% DMSO sensitivity']


# In[20]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)


# In[37]:


dataset_ids = [16459, 16460]
datasets = datasets.reindex(index=dataset_ids)


# In[29]:


original_data = original_data1['data'].to_frame().join(original_data2['data'], how='outer', lsuffix='_1', rsuffix='_2')


# In[31]:


original_data.columns = dataset_ids


# In[33]:


original_data[original_data.isnull()] = 0


# # Load & process tested strains

# In[35]:


tested = pd.read_excel('raw_data/DELETION LIBRARY.xlsx', sheet_name='DELETION LIBRARY', skiprows=1)
tested = np.array(tested['ORF name'].unique())
tested = clean_orf(tested)
tested = translate_sc(tested, to='orf')
tested[tested == 'YELOO1C'] = 'YEL001C'
# Make sure everything translated ok
t = looks_like_orf(tested)
print(tested[~np.array(t)])
tested = np.setdiff1d(tested, np.array(['YMR41W']))


# In[36]:


tested = pd.DataFrame(index=tested)


# # Prepare the final dataset

# In[38]:


data = tested.join(original_data, how='outer')


# In[40]:


data.columns = datasets.loc[data.columns,'name']


# In[44]:


data[data.isnull()] = 0


# In[45]:


data = data.groupby(data.index).mean()


# In[46]:


# Create row index
data.index.name='orf'


# In[47]:


print('Final data dimensions: %d x %d' % (data.shape))


# In[48]:


print((data<0).sum(axis=0))


# # Print out

# In[49]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[50]:


from IO.save_data_to_db2 import *


# In[51]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[52]:


save_data_to_db(data, paper_pmid)


# In[ ]:




