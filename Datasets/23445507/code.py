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


paper_pmid = 23445507
paper_name = 'troppens_morrissey_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/Table S1.xlsx', sheet_name='Sheet1')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data['orf'] = original_data['Systematic name'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[9]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[11]:


original_data.loc[original_data['orf']=='YDR369', 'orf'] = 'YDR369C'
original_data.loc[original_data['orf']=='YOR382', 'orf'] = 'YOR382W'


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[19]:


data_dict = {'+': 1,'++': 2,'-': -1, 'nd': np.nan,'+/-': 0, '(+)/-': 0}


# In[20]:


original_data['data'] = original_data['120 μg/ml DAPG']


# In[21]:


original_data['data'] = [data_dict[x] for x in original_data['data']]


# In[22]:


original_data.set_index('orf', inplace=True)


# # Prepare the final dataset

# In[25]:


dataset_ids = [16562]


# In[26]:


datasets = datasets.reindex(index=dataset_ids)


# In[27]:


data = original_data[['data']].copy()


# In[28]:


data.columns = datasets['name'].values


# In[29]:


data = data.groupby(data.index).mean()


# In[30]:


# Create row index
data.index.name='orf'


# In[31]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[32]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[33]:


from IO.save_data_to_db2 import *


# In[34]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[35]:


save_data_to_db(data, paper_pmid)


# In[ ]:




