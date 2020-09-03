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


paper_pmid = 31270132
paper_name = 'hoffert_strome_2019' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/Table_S1.xlsx', sheet_name='Sheet1', skiprows=2)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data['ORF name'] = original_data['ORF name'].astype(str)


# In[8]:


# Eliminate all white spaces & capitalize
original_data['ORF name'] = clean_orf(original_data['ORF name'])


# In[9]:


# Translate to ORFs 
original_data['ORF name'] = translate_sc(original_data['ORF name'], to='orf')


# In[11]:


typo_fixes = {'YCLO51W': 'YCL051W','YHR139C-': 'YHR139C-A','YGR122C-': 'YGR122C-A'}


# In[12]:


for orf in typo_fixes.keys():
    original_data.loc[original_data['ORF name']==orf,'ORF name'] = typo_fixes[orf]


# In[ ]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF name'])


# In[14]:


# Remove the 1's at the end of certain ORFs
for orf in original_data.loc[~t,'ORF name'].values:
    new_orf = orf.rstrip('1')
    original_data.loc[original_data['ORF name']==orf,'ORF name'] = new_orf


# In[15]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF name'])
print(original_data.loc[~t,])


# In[17]:


original_data = original_data.loc[t,]


# In[18]:


original_data.head()


# In[19]:


data_replacements = {0: 0, '+': 1, '++': 2, '+++': 3, '++++': 4}


# In[20]:


original_data = original_data[['ORF name','A','B','A.1','B.1','A.2','B.2','A.3','B.3']]


# In[21]:


original_data.set_index('ORF name', inplace=True)


# In[24]:


for c in original_data.columns:
    original_data[c+'_num'] = original_data[c].apply(lambda x: data_replacements[x] if x in data_replacements else np.nan)


# In[29]:


original_data['data'] = original_data[['A_num','B_num','A.1_num','B.1_num','A.2_num','B.2_num','A.3_num','B.3_num']].sum(axis=1)


# In[33]:


original_data.shape


# In[35]:


original_data['num_vals'] = original_data[['A_num','B_num','A.1_num','B.1_num','A.2_num','B.2_num','A.3_num','B.3_num']].apply(lambda x: ~np.isnan(x)).sum(axis=1)


# In[39]:


original_data = original_data.loc[original_data['num_vals']>0,]


# # Prepare the final dataset

# In[41]:


dataset_ids = [16548]


# In[42]:


datasets = datasets.reindex(index=dataset_ids)


# In[43]:


data = original_data['data'].to_frame()


# In[44]:


data.columns = datasets['name'].values


# In[45]:


data = data.groupby(data.index).mean()


# In[46]:


# Create row index
data.index.name='orf'


# In[47]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[51]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[52]:


from IO.save_data_to_db2 import *


# In[53]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[54]:


save_data_to_db(data, paper_pmid)


# In[ ]:




