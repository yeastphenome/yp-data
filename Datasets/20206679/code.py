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


paper_pmid = 20206679
paper_name = 'zhao_jiang_2010' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/hits.xlsx', sheet_name='Sheet1')


# In[6]:


original_data.columns = [x.strip() for x in original_data.columns]


# In[7]:


original_data.head()


# In[8]:


cols = ['Systemic Name','Systemic Name.1','Systemic Name.2','Systemic Name.3']
original_data = pd.concat([original_data[col] for col in cols], axis=0).to_frame()


# In[9]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[10]:


original_data[0] = original_data[0].astype(str)


# In[11]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data[0])


# In[12]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'].values, to='orf')


# In[13]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[14]:


original_data.drop(index=original_data.loc[~t,].index, inplace=True)


# In[15]:


original_data['data'] = -1


# In[16]:


original_data.set_index('orfs', inplace=True)


# # Load & process tested strains

# In[17]:


tested = pd.read_excel('raw_data/DELETION LIBRARY.xlsx', sheet_name='DELETION LIBRARY', skiprows=1)
tested = np.array(tested['ORF name'].unique())
tested = clean_orf(tested)
tested = translate_sc(tested, to='orf')
tested[tested == 'YELOO1C'] = 'YEL001C'
# Make sure everything translated ok
t = looks_like_orf(tested)
print(tested[~np.array(t)])
tested = np.setdiff1d(tested, np.array(['YMR41W']))


# In[18]:


tested = pd.DataFrame(index=tested)


# # Prepare the final dataset

# In[19]:


dataset_ids = [16456]


# In[20]:


datasets = datasets.reindex(index=dataset_ids)


# In[21]:


data = tested.join(original_data['data'], how='outer')


# In[22]:


data[data['data'].isnull()] = 0


# In[23]:


data = data.groupby(data.index).mean()


# In[24]:


# Create row index
data.index.name='orf'


# In[25]:


data.columns = datasets['name'].values


# In[26]:


print('Final data dimensions: %d x %d' % (data.shape))


# In[27]:


print((data<0).sum(axis=0))


# # Print out

# In[28]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[51]:


from IO.save_data_to_db2 import *


# In[52]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[53]:


save_data_to_db(data, paper_pmid)


# In[ ]:




