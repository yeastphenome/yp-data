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


paper_pmid = 18780730
paper_name = 'sinha_steinmetz_2008' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_csv('raw_data/deletion_pool_data.txt', sep='\t')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.columns


# In[8]:


original_data.head()


# In[9]:


original_data['orfs'] = original_data['orf::batch:tagtype'].apply(lambda x: x.split(':')[0])


# In[10]:


original_data['orfs'] = original_data['orfs'].astype(str)


# In[11]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[12]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[13]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[14]:


original_data['37C'] = (original_data['37C_T5'] / original_data['T0']) / (original_data['30C_T5'] / original_data['T0'])


# In[15]:


original_data.sort_values(by='37C', ascending=False)[['orfs','T0','30C_T1','30C_T2','30C_T3','30C_T4','30C_T5','37C_T1','37C_T2','37C_T3','37C_T4','37C_T5']].head()


# In[16]:


original_data['rapa_12h'] = original_data['30C_RAPA_T1'] / original_data['T0']
original_data['rapa_24h'] = original_data['30C_RAPA_T2'] / original_data['T0']
original_data['rapa_36h'] = original_data['30C_RAPA_T3'] / original_data['T0']


# In[17]:


original_data.set_index('orfs', inplace=True)


# In[18]:


# Splits homozygous and heterozygous mutants
original_data_hom = original_data.loc[original_data['zygosity']=='hom'].copy()
original_data_het = original_data.loc[original_data['zygosity']=='het'].copy()


# In[19]:


original_data_hom = original_data_hom.groupby(original_data_hom.index).mean()
original_data_hom.shape


# In[20]:


original_data_het = original_data_het.groupby(original_data_het.index).mean()
original_data_het.shape


# In[21]:


# Pull them back together
data = original_data_hom[['37C','rapa_12h','rapa_24h','rapa_36h']].join(original_data_het[['37C','rapa_12h','rapa_24h','rapa_36h']],
                                                                       lsuffix='_hom', rsuffix='_het', how='outer')


# # Prepare the final dataset

# In[22]:


dataset_ids = [16511,16639,16640,16641,16638,16642,16643,16644]


# In[23]:


datasets = datasets.reindex(index=dataset_ids)


# In[24]:


data.columns = datasets['name'].values


# In[25]:


data = data.groupby(data.index).mean()


# In[26]:


# Create row index
data.index.name='orf'


# In[27]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[28]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[50]:


from IO.save_data_to_db2 import *


# In[52]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[53]:


save_data_to_db(data, paper_pmid)


# In[ ]:




