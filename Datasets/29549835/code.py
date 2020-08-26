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


paper_pmid = 29549835
paper_name = 'ng_chia_2018' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[6]:


original_data = pd.read_excel('raw_data/Supplementary Worksheet (EJMC 2018).xlsx', sheet_name='HOP_peptide 50 ')


# In[7]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[8]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[10]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'], to='orf')


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[12]:


original_data = original_data.loc[t,:]


# In[13]:


original_data.set_index('ORF', inplace=True)


# # Prepare the final dataset

# In[14]:


dataset_ids = [16443]


# In[15]:


datasets = datasets.reindex(index=dataset_ids)


# In[16]:


data = original_data['logFC'].to_frame()


# In[17]:


data.columns = datasets['name'].values


# In[18]:


data = data.groupby(data.index).mean()


# In[19]:


# Create row index
data.index.name='orf'


# In[20]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[21]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[22]:


from IO.save_data_to_db2 import *


# In[23]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[24]:


save_data_to_db(data, paper_pmid)


# In[ ]:




