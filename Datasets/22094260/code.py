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


paper_pmid = 22094260
paper_name = 'skrtic_schimmer_2011' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_csv('raw_data/het_damp.rawsummary', sep='\t')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


# First, eliminate the data for the DAMP strains
original_data = original_data.loc[~original_data['Hybridization REF'].str.contains('DAMP'),]
original_data.shape


# In[9]:


# Now, extract the ORF
original_data['orf'] = original_data['Hybridization REF'].apply(lambda x: x[0:x.find(':')])


# In[10]:


original_data['orf'] = original_data['orf'].astype(str)


# In[11]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[12]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[13]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[14]:


original_data = original_data.loc[t,]


# In[15]:


original_data.set_index('orf', inplace=True)


# In[16]:


for c in original_data.columns.values[1:]:
    original_data[c] = pd.to_numeric(original_data[c])


# In[17]:


# Take the average of the 2 YPGE_DMSO controls
original_data['YPGE_DMSO_avg'] = original_data[['10_11_24_YPGE_DMSO','10_11_24_YPGE_DMSO_2']].mean(axis=1)


# In[18]:


# Divide each treatment by its control
original_data['10_11_24_YPGE_chloramph_0.79_norm'] = original_data['10_11_24_YPGE_chloramph_0.79'] / original_data['YPGE_DMSO_avg']
original_data['10_11_24_YPGE_chloramph_0.99_norm'] = original_data['10_11_24_YPGE_chloramph_0.99'] / original_data['YPGE_DMSO_avg']
original_data['10_11_24_YPGE_doxorub.12.5_norm'] = original_data['10_11_24_YPGE_doxorub.12.5'] / original_data['YPGE_DMSO_avg']
original_data['10_11_24_YPGE_linezol_47.1_norm'] = original_data['10_11_24_YPGE_linezol_47.1'] / original_data['YPGE_DMSO_avg']

original_data['10_12_10_tigecyc51.5uM_norm'] = original_data['10_12_10_tigecyc51.5uM'] / original_data['10_12_10_tigecycDMSOctrl']
original_data['10_12_10_tigecyc64.4uM_norm'] = original_data['10_12_10_tigecyc64.4uM'] / original_data['10_12_10_tigecycDMSOctrl']
original_data['10_12_10_tigecyc80.5uM_norm'] = original_data['10_12_10_tigecyc80.5uM'] / original_data['10_12_10_tigecycDMSOctrl']


# # Prepare the final dataset

# In[19]:


dataset_ids = [16572,16591,16570,16573,16571,16592,16593]


# In[20]:


datasets = datasets.reindex(index=dataset_ids)


# In[21]:


cols_to_keep = ['10_11_24_YPGE_chloramph_0.79_norm','10_11_24_YPGE_chloramph_0.99_norm',
                '10_11_24_YPGE_doxorub.12.5_norm',
                '10_11_24_YPGE_linezol_47.1_norm',
                '10_12_10_tigecyc51.5uM_norm','10_12_10_tigecyc64.4uM_norm','10_12_10_tigecyc80.5uM_norm']


# In[22]:


data = original_data.loc[:,cols_to_keep]


# In[24]:


data.columns = datasets['name'].values


# In[25]:


data = data.groupby(data.index).mean()


# In[26]:


# Create row index
data.index.name='orf'


# In[27]:


print('Final data dimensions: %d x %d' % (data.shape))


# In[28]:


data.head()


# # Print out

# In[39]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[40]:


from IO.save_data_to_db2 import *


# In[41]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[42]:


save_data_to_db(data, paper_pmid)


# In[ ]:




