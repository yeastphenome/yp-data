#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd

import sys

from os.path import expanduser
sys.path.append(expanduser('~') + '/Lab/Utils/Python/')

from Conversions.translate import *
from Strings.is_a import *


# # Initial setup

# In[3]:


paper_pmid = 20206679
paper_name = 'zhao_jiang_2010' 


# In[4]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[5]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[15]:


original_data = pd.read_excel('raw_data/hits.xlsx', sheet_name='Sheet1')


# In[16]:


original_data.columns = [x.strip() for x in original_data.columns]


# In[17]:


original_data.head()


# In[18]:


cols = ['Systemic Name','Systemic Name.1','Systemic Name.2','Systemic Name.3']
original_data = pd.concat([original_data[col] for col in cols], axis=0).to_frame()


# In[20]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[23]:


original_data[0] = original_data[0].astype(str)


# In[25]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data[0])


# In[26]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'].values, to='orf')


# In[28]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[29]:


original_data.drop(index=original_data.loc[~t,].index, inplace=True)


# In[30]:


original_data['data'] = -1


# In[31]:


original_data.set_index('orfs', inplace=True)


# # Load & process tested strains

# In[33]:


tested = pd.read_excel('raw_data/DELETION LIBRARY.xlsx', sheet_name='DELETION LIBRARY', skiprows=1)
tested = np.array(tested['ORF name'].unique())
tested = clean_orf(tested)
tested = translate_sc(tested, to='orf')
tested[tested == 'YELOO1C'] = 'YEL001C'
# Make sure everything translated ok
t = looks_like_orf(tested)
print(tested[~np.array(t)])
tested = np.setdiff1d(tested, np.array(['YMR41W']))


# In[37]:


tested = pd.DataFrame(index=tested)


# # Prepare the final dataset

# In[34]:


dataset_ids = [16456]


# In[35]:


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




