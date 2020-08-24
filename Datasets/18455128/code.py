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


paper_pmid = 18455128
paper_name = 'lain_westwood_2008' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/ranked list of hypersensitivity data.xlsx', sheet_name='Ranked list')


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


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[11]:


original_data['data'] = -1


# In[12]:


original_data.set_index('ORF', inplace=True)


# # Load & process tested strains

# In[15]:


tested = pd.read_excel('raw_data/strain list.xlsx', sheet_name='Strains')


# In[16]:


tested = tested['ORF']


# In[17]:


tested = clean_orf(tested)


# In[18]:


tested = translate_sc(tested, to='orf')


# In[21]:


# Make sure everything translated ok
t = looks_like_orf(tested)
print(tested[~np.array(t),])


# In[23]:


missing = [orf for orf in original_data.index.values if orf not in tested]
missing


# In[26]:


tested = np.append(tested, missing)


# In[28]:


tested = np.unique(tested)
tested.shape


# # Prepare the final dataset

# In[29]:


dataset_ids = [16497]


# In[30]:


datasets = datasets.reindex(index=dataset_ids)


# In[31]:


data = pd.DataFrame(index=tested, columns=datasets['name'].values, data=0)


# In[32]:


data.loc[original_data.index, datasets['name'].values[0]] = original_data['data']


# In[33]:


data = data.groupby(data.index).mean()


# In[34]:


# Create row index
data.index.name='orf'


# In[35]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[37]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[38]:


from IO.save_data_to_db2 import *


# In[39]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[40]:


save_data_to_db(data, paper_pmid)


# In[ ]:




