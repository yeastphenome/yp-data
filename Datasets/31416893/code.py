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


paper_pmid = 31416893
paper_name = 'bhat_sengupta_2019' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[23]:


original_data = pd.read_excel('raw_data/whole screen.xlsx', sheet_name='Sheet2')


# In[24]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[25]:


original_data['Name'] = original_data['Name'].astype(str)


# In[26]:


# Eliminate all white spaces & capitalize
original_data['Name'] = clean_genename(original_data['Name'])


# In[27]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['Name'], to='orf')


# In[28]:


name_fix_map = {'GON2': 'YLL033W','CRS5':'YOR031W','SRF5':'YOR041C','MOR1':'YDR366C','BOP1':'YPL221W','FMO': 'YHR176W','GON3':'YHR177W','FMP53':'YLR201C',
               'OCT':'YKL134C','FMP17':'YGR033C','FMP31':'YOR286W','SRF6':'YNL179C','SWS1':'YDR290W','AAD6':'YFL056C','YSN1':'YNR065C','FMP35':'YIL157C','SDL1':'YIL167W',
                'ZSP1':'YBR287W','SRF4':'YDL023C','FLO8':'YER109C'}


# In[29]:


for g in name_fix_map.keys():
    original_data.loc[original_data['orfs']==g,'orfs'] = name_fix_map[g]


# In[30]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[33]:


original_data['data'] = original_data['Cysteine'] / original_data['Control']


# In[34]:


original_data.set_index('orfs', inplace=True)


# # Prepare the final dataset

# In[35]:


dataset_ids = [16547]


# In[36]:


datasets = datasets.reindex(index=dataset_ids)


# In[37]:


data = original_data['ratio'].to_frame()


# In[38]:


data.columns = datasets['name'].values


# In[39]:


data = data.groupby(data.index).mean()


# In[40]:


# Create row index
data.index.name='orf'


# In[41]:


print('Final data dimensions: %d x %d' % (data.shape))


# In[46]:


data.sort_values(by = datasets['name'].values[0], ascending=True).head(n=20)


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




