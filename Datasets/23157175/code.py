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


paper_pmid = 23157175
paper_name = 'zhang_jiang_2013' 


# In[98]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[99]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[45]:


original_data1 = pd.read_excel('raw_data/hits.xlsx', sheet_name='DMSO4')
original_data2 = pd.read_excel('raw_data/hits.xlsx', sheet_name='DMSO8')


# In[46]:


original_data1.columns = [x.strip() for x in original_data1.columns]
original_data2.columns = [x.strip() for x in original_data2.columns]


# In[47]:


original_data2.head()


# In[48]:


original_data1.drop(index=original_data1.loc[original_data1['DMSO sensitivity'].isnull()].index, inplace=True)
original_data2.drop(index=original_data2.loc[original_data2['8% DMSO sensitivity'].isnull()].index, inplace=True)


# In[49]:


original_data1['DMSO sensitivity'] = original_data1['DMSO sensitivity'].apply(lambda x: len(x.strip()))
original_data2['8% DMSO sensitivity'] = original_data2['8% DMSO sensitivity'].apply(lambda x: len(x.strip()))


# In[50]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[51]:


original_data1['genes'] = original_data1['Gene'].astype(str)
original_data2['genes'] = original_data2['Gene'].astype(str)


# In[52]:


# Eliminate all white spaces & capitalize
original_data1['genes'] = clean_genename(original_data1['genes'])
original_data2['genes'] = clean_genename(original_data2['genes'])


# In[70]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['genes'].values, to='orf')
original_data2['orfs'] = translate_sc(original_data2['genes'].values, to='orf')


# In[71]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[72]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[73]:


original_data1['data'] = -original_data1['DMSO sensitivity']
original_data2['data'] = -original_data2['8% DMSO sensitivity']


# In[74]:


original_data1.set_index('orfs', inplace=True)
original_data2.set_index('orfs', inplace=True)


# # Load & process tested strains

# In[87]:


tested = pd.read_excel('raw_data/DELETION LIBRARY.xlsx', sheet_name='DELETION LIBRARY', skiprows=1)


# In[88]:


tested = np.array(tested['ORF name'].unique())


# In[89]:


tested = clean_orf(tested)


# In[90]:


tested = translate_sc(tested, to='orf')


# In[84]:


tested[tested == 'YELOO1C'] = 'YEL001C'


# In[85]:


# Make sure everything translated ok
t = looks_like_orf(tested)
print(tested[~np.array(t)])


# In[91]:


tested = np.setdiff1d(tested, np.array(['YMR41W']))


# # Prepare the final dataset

# In[100]:


dataset_ids = [16459, 16460]


# In[101]:


datasets = datasets.reindex(index=dataset_ids)


# In[102]:


data = pd.DataFrame(index=tested, columns=datasets['name'].values, data=0)


# In[103]:


data.loc[original_data1.index, datasets['name'].values[0]] = original_data1['data']
data.loc[original_data2.index, datasets['name'].values[1]] = original_data2['data']


# In[104]:


data = data.groupby(data.index).mean()


# In[105]:


# Create row index
data.index.name='orf'


# In[106]:


print('Final data dimensions: %d x %d' % (data.shape))


# In[108]:


print((data<0).sum(axis=0))


# # Print out

# In[109]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[110]:


from IO.save_data_to_db2 import *


# In[111]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[112]:


save_data_to_db(data, paper_pmid)


# In[ ]:




