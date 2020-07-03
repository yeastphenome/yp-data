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


paper_pmid = 19503795
paper_name = 'westmoreland_bennett_2009' 


# In[4]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[5]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[23]:


original_data1 = pd.read_excel('raw_data/journal.pone.0005830.s001.XLS', sheet_name='Sheet1')
original_data2 = pd.read_excel('raw_data/journal.pone.0005830.s002.XLS', sheet_name='Sheet1')


# In[24]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[25]:


original_data1['ORF'] = original_data1['ORF'].astype(str)
original_data2['ORF'] = original_data2['ORF'].astype(str)


# In[26]:


# Eliminate all white spaces & capitalize
original_data1['ORF'] = clean_orf(original_data1['ORF'])
original_data2['ORF'] = clean_orf(original_data2['ORF'])


# In[27]:


# Translate to ORFs 
original_data1['ORF'] = translate_sc(original_data1['ORF'], to='orf')
original_data2['ORF'] = translate_sc(original_data2['ORF'], to='orf')


# In[31]:


# Make sure everything translated ok
t1 = looks_like_orf(original_data1['ORF'])
print(original_data1.loc[~t1,])


# In[32]:


# Make sure everything translated ok
t2 = looks_like_orf(original_data2['ORF'])
print(original_data2.loc[~t2,])


# In[33]:


data1 = original_data1.loc[t1,['ORF','DoxS (2n)']]
data2 = original_data2.loc[t2,['ORF','DoxS (2n)']]


# In[34]:


# "3" denotes the most sensitive strains; in addition, shifting all values by 2 to make them more extreme and unify with data2 (less sensitive phenotypes)
data1['DoxS (2n)'] = -data1['DoxS (2n)']-2


# In[35]:


data2['DoxS (2n)'] = -data2['DoxS (2n)']


# In[40]:


data1.set_index('ORF', inplace=True)
data2.set_index('ORF', inplace=True)


# In[44]:


data = data1.join(data2, lsuffix='_1', rsuffix='_2', how='outer')


# In[47]:


data['data'] = data[['DoxS (2n)_1','DoxS (2n)_2']].mean(axis=1)


# In[51]:


data.drop(columns=['DoxS (2n)_1','DoxS (2n)_2'], inplace=True)


# # Prepare the final dataset

# In[52]:


dataset_ids = [16454]


# In[53]:


datasets = datasets.reindex(index=dataset_ids)


# In[54]:


data.columns = datasets['name'].values


# In[55]:


data = data.groupby(data.index).mean()


# In[56]:


# Create row index
data.index.name='orf'


# In[57]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[62]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[63]:


from IO.save_data_to_db2 import *


# In[64]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[65]:


save_data_to_db(data, paper_pmid)


# In[ ]:




