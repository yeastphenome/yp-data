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


paper_pmid = 32469861
paper_name = 'liu_liu_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:





# In[6]:


print('Original data dimensions: %d x %d' % (tested_strains.shape))


# In[7]:


original_data = pd.read_excel('raw_data/journal.pgen.1008798.s007.xlsx', sheet_name='Combined data', skiprows=2)


# In[8]:


original_data.columns = ['rank1','orf1','','rank2','orf2','','rank3','orf3','']


# In[9]:


original_data.head()


# In[10]:


ranks = pd.concat([original_data['rank1'], original_data['rank2'], original_data['rank3']], axis=0, ignore_index=True)


# In[11]:


orfs = pd.concat([original_data['orf1'], original_data['orf2'], original_data['orf3']], axis=0, ignore_index=True)


# In[12]:


hit_data = ranks.to_frame().join(orfs.to_frame(), how='outer', lsuffix='_rank', rsuffix='_orf')


# In[13]:


hit_data.head()


# In[14]:


hit_data['0_orf'] = hit_data['0_orf'].astype(str)


# In[15]:


# Eliminate all white spaces & capitalize
hit_data['0_orf'] = clean_orf(hit_data['0_orf'])


# In[17]:


# Translate to ORFs 
hit_data['orfs'] = translate_sc(hit_data['0_orf'], to='orf')


# In[18]:


# Make sure everything translated ok
t = looks_like_orf(hit_data['orfs'])
print(hit_data.loc[~t,])


# In[19]:


hit_data = hit_data.loc[t,:]


# In[20]:


hit_data.set_index('orfs', inplace=True)


# In[22]:


hit_data['data'] = 1


# # Load & process tested strains

# In[23]:


tested_strains = pd.read_excel('raw_data/Original data after SGA Scoring sorted.xls', sheet_name='Combined data')


# In[24]:


tested_strains['Array ORF'] = tested_strains['Array ORF'].astype(str)


# In[25]:


# Eliminate all white spaces & capitalize
tested_strains['Array ORF'] = clean_orf(tested_strains['Array ORF'])


# In[26]:


# Translate to ORFs 
tested_strains['orfs'] = translate_sc(tested_strains['Array ORF'], to='orf')


# In[27]:


# Make sure everything translated ok
t = looks_like_orf(tested_strains['orfs'])
print(tested_strains.loc[~t,])


# In[28]:


tested = tested_strains['orfs'].unique()


# In[36]:


missing = [orf for orf in hit_data.index.values if orf not in tested]


# In[37]:


missing


# # Prepare the final dataset

# In[29]:


dataset_ids = [16543]


# In[30]:


datasets = datasets.reindex(index=dataset_ids)


# In[31]:


data = pd.DataFrame(index=tested, columns=datasets['name'].values, data=0)


# In[32]:


data.loc[hit_data.index, datasets['name'].values[0]] = hit_data['data']


# In[33]:


data = data.groupby(data.index).mean()


# In[34]:


# Create row index
data.index.name='orf'


# In[35]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[38]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[41]:


from IO.save_data_to_db2 import *


# In[42]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[43]:


save_data_to_db(data, paper_pmid)


# In[ ]:




