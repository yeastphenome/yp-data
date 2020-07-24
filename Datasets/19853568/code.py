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


paper_pmid = 19853568
paper_name = 'carroll_drubin_2009' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[13]:


original_data = pd.read_excel('raw_data/1-s2.0-S1534580709003438-mmc2-2.xls', sheen_name='Table S1', skiprows=7, header=None)


# In[14]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[16]:


original_data = original_data.loc[original_data[1]=='KO']


# In[17]:


original_data.head()


# In[33]:


original_data['ORF'] = original_data[3].astype(str)


# In[34]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[35]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'].values, to='orf')


# In[37]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[40]:


original_data[4] = [s.strip() for s in original_data[4]]


# In[41]:


data_vals = original_data[4].unique()


# In[43]:


data_vals_dict = {'++': 2, '+': 1, '-': -1, '- -': -2}


# In[44]:


original_data['data'] = [data_vals_dict[v] for v in original_data[4]]


# In[46]:


original_data.set_index('ORF', inplace=True)


# # Prepare the final dataset

# In[47]:


dataset_ids = [16455]


# In[48]:


datasets = datasets.reindex(index=dataset_ids)


# In[50]:


data = original_data[['data']]


# In[52]:


data.columns = datasets['name'].values


# In[53]:


data = data.groupby(data.index).mean()


# In[54]:


# Create row index
data.index.name='orf'


# In[55]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[56]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[57]:


from IO.save_data_to_db2 import *


# In[58]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[59]:


save_data_to_db(data, paper_pmid)


# In[ ]:




