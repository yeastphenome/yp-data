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


paper_pmid = 19324962
paper_name = 'holbein_dichtl_2009' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/Holbein_Data_S2.xlsx', sheet_name='Feuil1')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[9]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'], to='orf')


# In[10]:


original_data = original_data.loc[original_data['med-12'].notnull() | original_data['med-18'].notnull(),:]


# In[11]:


original_data.shape


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[13]:


original_data.set_index('ORF', inplace=True)


# In[14]:


original_data = original_data[['med-12','med-18']]


# # Prepare the final dataset

# In[15]:


dataset_ids = [16452,16453]


# In[16]:


datasets = datasets.reindex(index=dataset_ids)


# In[17]:


data = original_data.copy()


# In[18]:


data.columns = datasets['name'].values


# In[19]:


data = data.groupby(data.index).mean()


# In[20]:


# Create row index
data.index.name='orf'


# In[21]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[23]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[24]:


from IO.save_data_to_db2 import *


# In[25]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[26]:


save_data_to_db(data, paper_pmid)


# In[ ]:



