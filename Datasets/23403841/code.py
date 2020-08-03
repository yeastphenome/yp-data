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


paper_pmid = 23403841
paper_name = 'oconnor_vulpe_2012' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[59]:


original_data = {}
time = {}
for dt in np.arange(4)+1:
    original_data[dt] = pd.read_excel('raw_data/data sheet ' + str(dt) + '.xlsx', sheet_name='Sheet1', skiprows=1)
    time[dt] = original_data[dt].loc[0,'Treatment']
    original_data[dt].columns = ['genename','orf','82.5 mM','165 mM','330 mM']
    original_data[dt] = original_data[dt].drop(index=[0,1,2])
    
    original_data[dt]['orf'] = original_data[dt]['orf'].astype(str)
    original_data[dt]['orf'] = clean_orf(original_data[dt]['orf'])
    original_data[dt]['orf'] = translate_sc(original_data[dt]['orf'], to='orf')
    
    t = looks_like_orf(original_data[dt]['orf'])
    print(original_data[dt].loc[~t,])
    
    original_data[dt].set_index('orf', inplace=True)


# In[82]:


for dt in np.arange(4)+1:
    t = original_data[dt].copy()
    t.drop(columns=['genename'], inplace=True)
    t.columns = [c+'_'+time[dt] for c in t.columns]
    if dt == 1:
        data_5g = t.copy()
    elif dt == 2:
        data_15g = t.copy()
    elif dt == 3:
        data_5g = pd.concat((data_5g, t), axis=0)
    elif dt == 4:
        data_15g = pd.concat((data_15g, t), axis=0)


# In[83]:


data_5g[data_5g.isnull()] = 0
data_15g[data_15g.isnull()] = 0


# In[84]:


data_5g = data_5g.astype(float)
data_15g = data_15g.astype(float)


# In[85]:


data_5g = data_5g.groupby(data_5g.index).mean()
data_15g = data_15g.groupby(data_15g.index).mean()


# In[93]:


data = data_5g.join(data_15g, how='outer')
data[data.isnull()] = 0


# # Prepare the final dataset

# In[97]:


dataset_ids = [16531,16529,16526,16530,16528,16527]


# In[98]:


datasets = datasets.reindex(index=dataset_ids)


# In[99]:


data.columns = datasets['name'].values


# In[100]:


# Create row index
data.index.name='orf'


# In[101]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[103]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[104]:


from IO.save_data_to_db2 import *


# In[105]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[106]:


save_data_to_db(data, paper_pmid)


# In[ ]:




