#!/usr/bin/env python
# coding: utf-8

# In[2]:


import numpy as np
import pandas as pd

import sys

from os.path import expanduser
sys.path.append(expanduser('~') + '/Lab/Utils/Python/')

from Conversions.translate import *
from Strings.is_a import *


# # Initial setup

# In[34]:


paper_pmid = 18780730
paper_name = 'sinha_steinmetz_2008' 


# In[35]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[36]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[3]:


original_data = pd.read_csv('raw_data/deletion_pool_data.txt', sep='\t')


# In[4]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[6]:


original_data.columns


# In[7]:


original_data.head()


# In[10]:


original_data['orfs'] = original_data['orf::batch:tagtype'].apply(lambda x: x.split(':')[0])


# In[12]:


original_data['orfs'] = original_data['orfs'].astype(str)


# In[13]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[14]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[15]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[16]:


original_data['37C'] = (original_data['37C_T5'] / original_data['T0']) / (original_data['30C_T5'] / original_data['T0'])


# In[22]:


original_data.sort_values(by='37C', ascending=False)[['orfs','T0','30C_T1','30C_T2','30C_T3','30C_T4','30C_T5','37C_T1','37C_T2','37C_T3','37C_T4','37C_T5']].head()


# In[23]:


original_data['rapa_12h'] = original_data['30C_RAPA_T1'] / original_data['T0']
original_data['rapa_24h'] = original_data['30C_RAPA_T2'] / original_data['T0']
original_data['rapa_36h'] = original_data['30C_RAPA_T3'] / original_data['T0']


# In[24]:


original_data.set_index('orfs', inplace=True)


# In[25]:


# Splits homozygous and heterozygous mutants
original_data_hom = original_data.loc[original_data['zygosity']=='hom'].copy()
original_data_het = original_data.loc[original_data['zygosity']=='het'].copy()


# In[28]:


original_data_hom = original_data_hom.groupby(original_data_hom.index).mean()
original_data_hom.shape


# In[29]:


original_data_het = original_data_het.groupby(original_data_het.index).mean()
original_data_het.shape


# In[42]:


# Pull them back together
data = original_data_hom[['37C','rapa_12h','rapa_24h','rapa_36h']].join(original_data_het[['37C','rapa_12h','rapa_24h','rapa_36h']],
                                                                       lsuffix='_hom', rsuffix='_het', how='outer')


# # Prepare the final dataset

# In[43]:


dataset_ids = [16511,16639,16640,16641,16638,16642,16643,16644]


# In[44]:


datasets = datasets.reindex(index=dataset_ids)


# In[45]:


data.columns = datasets['name'].values


# In[46]:


data = data.groupby(data.index).mean()


# In[47]:


# Create row index
data.index.name='orf'


# In[48]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[51]:


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




