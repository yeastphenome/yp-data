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


paper_pmid = 26068574
paper_name = 'li_breeden_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[6]:


original_data = pd.read_excel('raw_data/table_1.xlsx', sheet_name='Sheet1')


# In[7]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[9]:


original_data['orfs'] = original_data['ORF ID'].astype(str)


# In[10]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[11]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[13]:


original_data = original_data.loc[t,]


# In[14]:


original_data['data'] = -1


# In[15]:


original_data.set_index('orfs', inplace=True)


# # Load & process tested strains

# In[21]:


tested = pd.read_excel('raw_data/library screen.xlsx', sheet_name='Sheet1')


# In[22]:


tested['orfs'] = tested['ORF name'].astype(str)
tested['orfs'] = clean_orf(tested['orfs'])


# In[23]:


tested.loc[tested['orfs']=='YLR287-A','orfs'] = 'YLR287C-A'


# In[24]:


tested['orfs'] = translate_sc(tested['orfs'], to='orf')


# In[25]:


# Make sure everything translated ok
t = looks_like_orf(tested['orfs'])
print(tested.loc[~t,])


# In[26]:


tested = tested.loc[t,]


# In[27]:


tested_orfs = tested['orfs'].unique()


# In[32]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]


# In[34]:


print(missing)


# In[35]:


# Adding the 2 missing hits to the list of tested strains
tested_orfs = np.array(list(tested_orfs) + missing)


# # Prepare the final dataset

# In[36]:


dataset_ids = [16601]


# In[37]:


datasets = datasets.reindex(index=dataset_ids)


# In[38]:


data = pd.DataFrame(index=tested_orfs, columns=datasets['name'].values, data=0)


# In[39]:


data.loc[original_data.index, datasets['name'].values[0]] = original_data['data']


# In[40]:


data = data.groupby(data.index).mean()


# In[41]:


# Create row index
data.index.name='orf'


# In[42]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[44]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[45]:


from IO.save_data_to_db2 import *


# In[46]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[47]:


save_data_to_db(data, paper_pmid)


# In[ ]:




