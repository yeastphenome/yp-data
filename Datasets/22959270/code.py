#!/usr/bin/env python
# coding: utf-8

# In[53]:


import numpy as np
import pandas as pd

from itertools import compress

import sys

from os.path import expanduser
sys.path.append(expanduser('~') + '/Lab/Utils/Python/')

from Conversions.translate import *
from Strings.is_a import *


# # Initial setup

# In[2]:


paper_pmid = 22959270
paper_name = 'zhang_lobachev_2012' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, names=['genes'])


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data['genes'] = original_data['genes'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[9]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[10]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[11]:


original_data['data'] = 1


# In[12]:


original_data.set_index('orfs', inplace=True)


# # Load & process tested strains

# In[55]:


tested = pd.read_excel('raw_data/mat_a_101501.xlsx', sheet_name='mat_a_101501', skiprows=1)


# In[56]:


tested = tested['ORF name'].unique()


# In[57]:


tested = tested.astype(str)


# In[58]:


tested = clean_orf(tested)


# In[59]:


tested = translate_sc(tested, to='orf')


# In[60]:


tested[tested == ['YLR287-A']] = 'YLR287C-A'


# In[61]:


tested = list(compress(tested, ~(tested=='NAN')))


# In[62]:


# Make sure everything translated ok
t = looks_like_orf(tested)


# In[63]:


list(compress(tested, ~t))


# In[65]:


# Test if all hits are present in tested
missing = [orf for orf in original_data.index.values if orf not in tested]
print(missing)


# # Prepare the final dataset

# In[66]:


dataset_ids = [16402]


# In[67]:


datasets = datasets.reindex(index=dataset_ids)


# In[68]:


data = pd.DataFrame(index=tested, columns=datasets['name'].values, data=0)


# In[69]:


data.loc[original_data.index, datasets['name'].values[0]] = original_data['data']


# In[70]:


data = data.groupby(data.index).mean()


# In[71]:


# Create row index
data.index.name='orf'


# In[72]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[75]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[76]:


from IO.save_data_to_db2 import *


# In[77]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[78]:


save_data_to_db(data, paper_pmid)


# In[ ]:




