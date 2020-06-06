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


# In[2]:


get_ipython().system('pwd')


# # Initial setup

# In[3]:


paper_pmid = 25773006
paper_name = 'du_jiang_2015' 


# In[4]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[5]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[26]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, names=['genes'])


# In[27]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[28]:


original_data['genes'] = original_data['genes'].astype(str)


# In[29]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[30]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[31]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])


# In[32]:


print(original_data.loc[~t,])


# In[33]:


original_data['data'] = -1


# In[34]:


original_data.set_index('orfs', inplace=True)


# In[35]:


dataset_ids = [16461]
datasets = datasets.reindex(index=dataset_ids)


# In[36]:


original_data = original_data['data'].to_frame()
original_data.columns = dataset_ids


# # Load & process tested strains

# In[15]:


tested = pd.read_excel('raw_data/DELETION LIBRARY.xlsx', sheet_name='DELETION LIBRARY', skiprows=1)
tested = np.array(tested['ORF name'].unique())
tested = clean_orf(tested)
tested = translate_sc(tested, to='orf')
tested[tested == 'YELOO1C'] = 'YEL001C'
# Make sure everything translated ok
t = looks_like_orf(tested)
print(tested[~np.array(t)])
tested = np.setdiff1d(tested, np.array(['YMR41W']))


# In[16]:


tested = pd.DataFrame(index=tested)


# # Prepare the final dataset

# In[38]:


data = tested.join(original_data, how='outer')


# In[39]:


data.columns = datasets.loc[data.columns,'name']


# In[19]:


data = data.groupby(data.index).mean()


# In[20]:


# Create row index
data.index.name='orf'


# In[21]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[22]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[23]:


from IO.save_data_to_db2 import *


# In[24]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[25]:


save_data_to_db(data, paper_pmid)


# In[ ]:




