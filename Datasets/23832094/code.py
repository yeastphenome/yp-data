#!/usr/bin/env python
# coding: utf-8

# In[8]:


import numpy as np
import pandas as pd

import sys

from os.path import expanduser
sys.path.append(expanduser('~') + '/Lab/Utils/Python/')

from Conversions.translate import *
from Strings.is_a import *


# # Initial setup

# In[9]:


paper_pmid = 23832094
paper_name = 'tun_wu_2013' 


# In[10]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[11]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[29]:


original_data = pd.read_excel('raw_data/c3mt00083d.xlsx', sheet_name='Sheet1', skiprows=1)


# In[30]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[31]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[32]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[33]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'], to='orf')


# In[34]:


original_data.loc[original_data['ORF']=='JR078W','ORF'] = 'YJR078W'
original_data.loc[original_data['ORF']=='YOR205CHOMDIP','ORF'] = 'YOR205C'


# In[35]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[36]:


original_data.set_index('ORF', inplace=True)


# In[37]:


# Normalize to Control
original_data.iloc[:,1:] = original_data.iloc[:,1:].div(original_data.iloc[:,0], axis=0)


# In[38]:


# Normalize to WT
original_data = original_data.div(original_data.iloc[0,:], axis=1)


# In[40]:


# Remove the WT
original_data = original_data.drop(index='BY4743')


# # Prepare the final dataset

# In[43]:


dataset_ids = [16483, 16481, 16479, 16482]


# In[44]:


datasets = datasets.reindex(index=dataset_ids)


# In[45]:


data = original_data.copy()


# In[46]:


data.columns = datasets['name'].values


# In[47]:


data = data.groupby(data.index).mean()


# In[48]:


# Create row index
data.index.name='orf'


# In[49]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[51]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[52]:


from IO.save_data_to_db2 import *


# In[53]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[54]:


save_data_to_db(data, paper_pmid)


# In[ ]:




