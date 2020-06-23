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

# In[61]:


paper_pmid = 24926745
paper_name = 'tun_wu_2014' 


# In[62]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[63]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[8]:


original_data = pd.read_excel('raw_data/c4mt00116h1.xlsx', sheet_name='Sheet1', skiprows=1)


# In[9]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[10]:


original_data.head()


# In[11]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[12]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[13]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'], to='orf')


# In[15]:


original_data.loc[original_data['ORF'].str.startswith('YOR205CHOMDIP'),'ORF'] = 'YOR205C'


# In[16]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[43]:


data = original_data[['ORF','Control','Al 1.6 mM','Al  3.2 mM']].copy()


# In[44]:


data.set_index('ORF', inplace=True)


# In[45]:


data['Control'] = pd.to_numeric(data['Control'], errors='coerce')
data['Al 1.6 mM'] = pd.to_numeric(data['Al 1.6 mM'], errors='coerce')
data['Al  3.2 mM'] = pd.to_numeric(data['Al  3.2 mM'], errors='coerce')


# In[48]:


data = data.div(data.loc['BY4743',:])


# In[49]:


data['Al 1.6 mM'] = data['Al 1.6 mM'] / data['Control']


# In[50]:


data['Al  3.2 mM'] = data['Al  3.2 mM'] / data['Control']


# In[55]:


data.drop(index='BY4743', inplace=True)


# In[54]:


data.sort_values(by='Al  3.2 mM', ascending=False).head()


# In[56]:


data = data.groupby(data.index).mean()


# In[57]:


data.shape


# # Prepare the final dataset

# In[58]:


dataset_ids = [16509,16477,16478]


# In[64]:


datasets = datasets.reindex(index=dataset_ids)


# In[66]:


data.columns = datasets['name'].values


# In[67]:


# Create row index
data.index.name='orf'


# In[68]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[69]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[70]:


from IO.save_data_to_db2 import *


# In[71]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[72]:


save_data_to_db(data, paper_pmid)


# In[ ]:




