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


paper_pmid = 20429770
paper_name = 'mir_rashed_smith_2010' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/Ech sup Table work 5.xls', sheet_name='Sheet1', skiprows=4)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[10]:


original_data.drop(columns=original_data.columns[2::3], inplace=True)


# In[17]:


cols_orfs = np.arange(32)[0::2]


# In[22]:


all_orfs = []
all_datasets = []
for ix_dataset in np.arange(16):
    col_orf = cols_orfs[ix_dataset]
    all_orfs = all_orfs + original_data.iloc[:,col_orf].to_list()
    all_datasets = all_datasets + [original_data.columns.values[col_orf]]


# In[23]:


all_orfs = np.unique(np.array(all_orfs))


# In[35]:


data = pd.DataFrame(index=all_orfs, columns=all_datasets)


# In[36]:


for ix_dataset in np.arange(16):
    col_orf = cols_orfs[ix_dataset]
    col_data = col_orf+1
    
    this_dataset_name = original_data.columns.values[col_orf]
    this_orfs = original_data.iloc[:, col_orf].astype(str).values[1:]
    this_data = original_data.iloc[:, col_data].values[1:]
    
    ix = (this_orfs != 'nan')
    this_data = this_data[ix]
    this_orfs = this_orfs[ix]
    
    data.loc[this_orfs, this_dataset_name] = this_data


# In[38]:


orfs = data.index.values


# In[39]:


# Eliminate all white spaces & capitalize
orfs = clean_orf(orfs)


# In[40]:


# Translate to ORFs 
orfs = translate_sc(orfs, to='orf')


# In[42]:


# Make sure everything translated ok
t = looks_like_orf(orfs)
print(orfs[~np.array(t),])


# In[43]:


orfs = orfs[np.array(t)]


# In[44]:


data = data.reindex(index=orfs)


# In[47]:


# Transform data to follow convention (0 = WT, negative values = less growth than WT)
data = -data
data[data.isnull()] = 0


# In[60]:


data = data.astype(float)


# In[48]:


data.head()


# # Prepare the final dataset

# In[51]:


dataset_names = dict.fromkeys(data.columns.values)


# In[53]:


dataset_names['SG1.Dark'] = 16457
dataset_names['SG1.UV'] = 16575
dataset_names['SGEPR.Dark'] = 16576
dataset_names['SGEPR.UV'] = 16577
dataset_names['SGEPF.Dark'] = 16578
dataset_names['SGEPF.UV'] = 16579
dataset_names['SGEPH.Dark'] = 16580
dataset_names['SGEPH.UV'] = 16581
dataset_names['SGEPLS.Dark'] = 16582
dataset_names['SGEPLS.UV'] = 16583
dataset_names['SG7a.Dark'] = 16584
dataset_names['SG7a.UV'] = 16585
dataset_names['SGEARa.Dark'] = 16586
dataset_names['SGEARa.UV'] = 16587
dataset_names['SGEARb.Dark'] = 16588
dataset_names['SGEARb.UV'] = 16589


# In[55]:


dataset_ids = [dataset_names[d] for d in data.columns.values]


# In[57]:


datasets = datasets.reindex(index=dataset_ids)


# In[58]:


data.columns = datasets['name'].values


# In[61]:


data = data.groupby(data.index).mean()


# In[62]:


# Create row index
data.index.name='orf'


# In[63]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[65]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[66]:


from IO.save_data_to_db2 import *


# In[67]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[68]:


save_data_to_db(data, paper_pmid)


# In[ ]:




