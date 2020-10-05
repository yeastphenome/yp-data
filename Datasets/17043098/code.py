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

# In[3]:


paper_pmid = 17043098
paper_name = 'doostzadeh_langston_2007' 


# In[4]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[5]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[8]:


original_data = pd.read_excel('raw_data/kfl131supp.xls', sheet_name='ToxSci Supplementary Data files', skiprows=1)


# In[9]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[10]:


original_data['strain'] = original_data['strain'].astype(str)


# In[11]:


mpp_columns = ['z_result_nq:04_10_28_19:mpp+:250:ug/ml::::20:hom_09_02',
              'z_result_nq:04_10_28_25:mpp+:250:ug/ml::::20:hom_09_02',
              'z_result_nq:04_11_04_06:mpp+:250:ug/ml::::20:hom_09_02']

paraquat_columns = ['z_result_nq:04_10_28_21:paraquat:5000:uM::::20:hom_09_02',
                   'z_result_nq:04_10_28_27:paraquat:5000:uM::::20:hom_09_02',
                   'z_result_nq:04_11_04_08:paraquat:5000:uM::::20:hom_09_02']


# In[12]:


# Extract ORF from string
orfs = original_data['strain'].apply(lambda x: x.split(':')[0])


# In[14]:


original_data['orfs'] = orfs


# In[15]:


# Eliminate all white spaces & capitalize
original_data['orfs'] = clean_orf(original_data['orfs'])


# In[16]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['orfs'], to='orf')


# In[17]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[18]:


original_data['mpp'] = original_data[mpp_columns].mean(axis=1)


# In[19]:


original_data['paraquat'] = original_data[paraquat_columns].mean(axis=1)


# In[20]:


original_data.set_index('orfs', inplace=True)


# # Prepare the final dataset

# In[37]:


dataset_ids = [16616,16615]


# In[38]:


datasets = datasets.reindex(index=dataset_ids)


# In[39]:


data = original_data[['mpp','paraquat']].copy()


# In[40]:


data.columns = datasets['name'].values


# In[41]:


data = data.groupby(data.index).mean()


# In[42]:


# The original fitness scores were calculated as UNT/TRT, so that high value = growth defect. 
# Our convention is the opposite (i.e., low value = growth defect), so need to flip the sign of the data.
data = -data


# In[43]:


# Create row index
data.index.name='orf'


# In[44]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[47]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[48]:


from IO.save_data_to_db2 import *


# In[49]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[50]:


save_data_to_db(data, paper_pmid)


# In[ ]:




