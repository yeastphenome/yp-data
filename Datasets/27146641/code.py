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


paper_pmid = 27146641
paper_name = 'johnson_wu_2016' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[62]:


original_data = pd.read_excel('raw_data/c6mt00039h1.xlsx', sheet_name='Sheet1', skiprows=1, 
                              names=['orf','h2o_t0','h2o_t16','chr5_t0','chr5_t16','chr1_t0','chr1_t16'])


# In[63]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[64]:


original_data.head()


# In[65]:


original_data['orf'] = original_data['orf'].astype(str)


# In[66]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[67]:


# If possible, fix typos, omissions, etc.
original_data.loc[original_data['orf'].str.contains('BY4743AVERAGEN128'),'orf'] = 'WT'


# In[68]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[69]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])


# In[70]:


print(original_data.loc[~t,])


# In[71]:


# Normalize by t16 by t0, treated vs untreated and mut vs wt
original_data['h2o_ratio'] = original_data['h2o_t16'] / original_data['h2o_t0']
original_data['chr5_ratio'] = original_data['chr5_t16'] / original_data['chr5_t0']
original_data['chr1_ratio'] = original_data['chr1_t16'] / original_data['chr1_t0']


# In[72]:


original_data['h2o_ratio_wt'] = original_data['h2o_ratio'] / original_data.loc[original_data['orf']=='WT','h2o_ratio'].values
original_data['chr5_ratio_wt'] = original_data['chr5_ratio'] / original_data.loc[original_data['orf']=='WT','chr5_ratio'].values
original_data['chr1_ratio_wt'] = original_data['chr1_ratio'] / original_data.loc[original_data['orf']=='WT','chr1_ratio'].values


# In[74]:


original_data['chr5_ratio_wt_unt'] = original_data['chr5_ratio_wt'] / original_data['h2o_ratio_wt']
original_data['chr1_ratio_wt_unt'] = original_data['chr1_ratio_wt'] / original_data['h2o_ratio_wt']


# In[75]:


original_data.head()


# In[79]:


# If the same strain is present more than once, average its values
data = original_data.groupby('orf')[['chr5_ratio_wt_unt','chr1_ratio_wt_unt']].mean()


# In[82]:


data.drop(index='WT', inplace=True)


# In[83]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Prepare the final dataset

# In[87]:


dataset_ids = [16447, 16446]


# In[88]:


datasets = datasets.reindex(index=dataset_ids)


# In[89]:


data.columns = datasets['name']


# In[90]:


# Create row index
data.index.name='orf'


# In[ ]:


data.


# # Print out

# In[20]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[21]:


from IO.save_data_to_db2 import *


# In[ ]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[22]:


save_data_to_db(data, paper_pmid)

