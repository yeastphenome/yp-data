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


paper_pmid = 26267134
paper_name = 'costa_texeira_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[14]:


original_data = pd.read_excel('raw_data/journal.pone.0135110.s002.XLSX',sheet_name='Resistance determinants', skiprows=4, header=None)


# In[15]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[16]:


original_data.head()


# In[17]:


original_data['orfs'] = original_data[1].astype(str)


# In[18]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[19]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[20]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[21]:


original_data = original_data.loc[t,]


# In[49]:


original_data['data'] = 1


# In[50]:


original_data.set_index('orfs', inplace=True)


# In[59]:


original_data.shape


# # Load & process tested strains

# In[24]:


tested1 = pd.read_excel('raw_data/BY4741-1stDelivery.xls', sheet_name='Tabelle1')
tested2 = pd.read_excel('raw_data/BY4741-2nd Delivery.xls', sheet_name='chr11_1yes')
tested3 = pd.read_excel('raw_data/BY4741-3rd Delivery.xls', sheet_name='Tabelle1')


# In[37]:


tested = pd.concat((tested1['ORF'], tested2['ORF'], tested3['orf']), axis=0)


# In[38]:


tested = clean_orf(tested)


# In[39]:


tested = translate_sc(tested, to='orf')


# In[40]:


tested = np.unique(tested)


# In[42]:


missing = [orf for orf in original_data['orfs'] if orf not in tested]


# In[44]:


tested = list(tested) + missing


# # Prepare the final dataset

# In[46]:


dataset_ids = [16462]


# In[47]:


datasets = datasets.reindex(index=dataset_ids)


# In[48]:


data = pd.DataFrame(index=tested, columns=datasets['name'].values, data=0)


# In[53]:


data.loc[original_data.index, datasets['name'].values[0]] = original_data['data']


# In[54]:


data = data.groupby(data.index).mean()


# In[55]:


# Create row index
data.index.name='orf'


# In[56]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[60]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[61]:


from IO.save_data_to_db2 import *


# In[62]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[63]:


save_data_to_db(data, paper_pmid)


# In[ ]:



