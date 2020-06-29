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

# In[3]:


paper_pmid = 24262822
paper_name = 'mattiazziusaj_petrovic_2014' 


# In[4]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[5]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[9]:


original_data = pd.read_csv('raw_data/MattiazziUsaj_Chemosphere_2014.csv', sep=',')


# In[10]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[11]:


original_data.head()


# In[12]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[13]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[14]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'], to='orf')


# In[15]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[26]:


data = original_data.loc[:,['ORF','Thiamethoxam','Actara','Acetamiprid','Mospilan','DMSO+Pyrrolidone','Confidor ']]


# In[27]:


data.set_index('ORF',inplace=True)


# # Prepare the final dataset

# In[28]:


dataset_ids = [16502,16506,16503,16505,16504,16507]


# In[29]:


datasets = datasets.reindex(index=dataset_ids)


# In[30]:


data.columns = datasets['name'].values


# In[31]:


data = data.groupby(data.index).mean()


# In[32]:


# Create row index
data.index.name='orf'


# In[33]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[35]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[36]:


from IO.save_data_to_db2 import *


# In[37]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[38]:


save_data_to_db(data, paper_pmid)


# In[ ]:




