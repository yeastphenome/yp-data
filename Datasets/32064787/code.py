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


paper_pmid = 32064787
paper_name = 'mattiazzi_usaj_andrews_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/msb199243-sup-0003-tableev2.xlsx', sheet_name='penetrance_phenotype_data')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[9]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])


# In[10]:


print(original_data.loc[~t,])


# In[11]:


# Eliminate strains that are not deletions
dels = original_data['StrainID'].str.startswith('DMA')
original_data = original_data.loc[dels.values,:]
print(original_data.shape)


# In[12]:


original_data.set_index('ORF', inplace=True)


# In[13]:


dataset_map = {'actin_aggregate': 16403,
               'actin_bright_patches': 16415,
               'actin_decreased_patch_number': 16416,
               'actin_depolarized_patches': 16417,
               'coat_aggregate': 16418,
               'coat_decreased_patch_number': 16419,
               'coat_depolarized_patches': 16420,
               'coat_increased_patch_number': 16421,
               'LE_fragmented': 16422,
               'LE_membrane': 16423,
               'LE_condensed': 16424,
               'vacuole_class_E': 16425,
               'vacuole_enlarged': 16426,
               'vacuole_VATPase_defect': 16427,
               'vacuole_fragmented': 16428,
               'vacuole_class_G': 16429,
               'vacuole_multilobed': 16430}


# In[14]:


# Get the relevant columns
original_data = original_data.loc[:, dataset_map.keys()]


# In[15]:


print(original_data.shape)


# In[25]:


# If the same strain is present more than once, average its values
data = original_data.groupby(original_data.index.values).mean()


# In[18]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Prepare the final dataset

# In[28]:


dataset_ids = [dataset_map[c] for c in data.columns.values]


# In[29]:


datasets = datasets.reindex(index=dataset_ids)


# In[ ]:


# Create row index
data.index.name='orf'
data.columns = datasets['name']


# # Print out

# In[35]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[36]:


from IO.save_data_to_db2 import *


# In[38]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[39]:


save_data_to_db(data, paper_pmid)


# In[ ]:




