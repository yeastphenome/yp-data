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


paper_pmid = 31638329
paper_name = 'valero_gonzalez_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[6]:


original_data = pd.read_excel('raw_data/LogFC todos.xlsx', sheet_name='Ctrol vs sulfito')


# In[7]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[9]:


orf_col = 'Systematic name'


# In[10]:


original_data[orf_col] = original_data[orf_col].astype(str)


# In[11]:


# Eliminate all white spaces & capitalize
original_data[orf_col] = clean_orf(original_data[orf_col])


# In[18]:


# Remove the trailing replicate numbers ("-1" and "-2")
def remove_trailing_rep(orf):
    reps = ['-1','-2']
    for r in reps:
        f = orf.find(r)
        if f > 0:
            orf = orf[:f]
    return orf


# In[21]:


original_data['orf'] = original_data[orf_col].apply(remove_trailing_rep)


# In[22]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[23]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[24]:


original_data['data'] = original_data['LogFC']


# In[25]:


original_data.set_index('orf', inplace=True)


# # Prepare the final dataset

# In[26]:


dataset_ids = [16544]


# In[27]:


datasets = datasets.reindex(index=dataset_ids)


# In[28]:


data = original_data[['data']].copy()


# In[29]:


data.columns = datasets['name'].values


# In[30]:


data = data.groupby(data.index).mean()


# In[31]:


# Create row index
data.index.name='orf'


# In[32]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[33]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[35]:


from IO.save_data_to_db2 import *


# In[36]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[37]:


save_data_to_db(data, paper_pmid)


# In[ ]:




