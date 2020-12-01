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


paper_pmid = 23378416
paper_name = 'zeidler_denfert_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[44]:


original_data = pd.read_excel('raw_data/Table1.xlsx', sheet_name='Sheet1')


# In[45]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[46]:


original_data['orfs'] = original_data['ORF'].astype(str)


# In[47]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[48]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[49]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[50]:


original_data = original_data.loc[t,:]


# In[51]:


original_data.set_index('orfs', inplace=True)


# # Prepare the final dataset

# In[52]:


dataset_ids = [16521,16519,16522]


# In[53]:


datasets = datasets.reindex(index=dataset_ids)


# In[54]:


data = original_data[['Colistin','Aminocandin','Combined']].copy()


# In[56]:


data['Colistin'] = pd.to_numeric(data['Colistin'], errors='coerce')
data['Aminocandin'] = pd.to_numeric(data['Aminocandin'], errors='coerce')
data['Combined'] = pd.to_numeric(data['Combined'], errors='coerce')


# In[57]:


data.columns = datasets['name'].values


# In[58]:


data = data.groupby(data.index).mean()


# In[59]:


# Create row index
data.index.name='orf'


# In[60]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[61]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[62]:


from IO.save_data_to_db2 import *


# In[63]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[64]:


save_data_to_db(data, paper_pmid)


# In[ ]:




