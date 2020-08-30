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


paper_pmid = 31422915
paper_name = 'barbosa_siniossoglou_2019' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/mmc2.xlsx', sheet_name='Sheet1', skiprows=7)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[8]:


d1 = 'ER+puncta'
d2 = 'Low signal'


# In[51]:


d1_ix = 0
d2_ix = original_data.loc[original_data.iloc[:,0]==d2].index.values[0]


# In[52]:


data1 = original_data.iloc[d1_ix:d2_ix,:].copy()
data2 = original_data.iloc[d2_ix+1:,:].copy()


# In[53]:


data1.columns = data1.iloc[0,:].copy()


# In[54]:


data2.columns = data2.iloc[0,:].copy()


# In[55]:


data1.head()


# In[56]:


c = 'Systematic Name'
data1[c] = data1[c].astype(str)
data2[c] = data2[c].astype(str)


# In[57]:


# Eliminate all white spaces & capitalize
data1[c] = clean_orf(data1[c])
data2[c] = clean_orf(data2[c])


# In[58]:


# Translate to ORFs 
data1[c] = translate_sc(data1[c], to='orf')
data2[c] = translate_sc(data2[c], to='orf')


# In[59]:


# Make sure everything translated ok
t = looks_like_orf(data1[c])
print(data1.loc[~t,])


# In[60]:


data1 = data1.loc[t,:]


# In[61]:


# Make sure everything translated ok
t = looks_like_orf(data2[c])
print(data2.loc[~t,])


# In[62]:


data2 = data2.loc[t,:]


# In[63]:


data1['data'] = 1
data2['data'] = 1


# In[64]:


data1.head()


# In[65]:


data1.set_index(c, inplace=True)
data2.set_index(c, inplace=True)


# In[96]:


data2.head()


# In[97]:


data = data1[['data']].join(data2[['data']], lsuffix='_1', rsuffix='_2', how='outer')


# In[98]:


data.loc[data['data_1'].isnull(),'data_1'] = 0
data.loc[data['data_2'].isnull(),'data_2'] = 0


# In[99]:


data.shape


# In[100]:


data.sum(axis=0)


# In[103]:


data2 = data.copy()


# # Load & process tested strains

# In[70]:


tested = pd.read_excel('raw_data/KO_DAmP_ORFs.xlsx', sheet_name='Sheet1', skiprows=1)


# In[71]:


tested = tested.iloc[:,0].to_frame()


# In[76]:


tested.columns = ['ORF']


# In[77]:


tested['ORF'] = clean_orf(tested['ORF'])


# In[78]:


tested['ORF'] = translate_sc(tested['ORF'], to='orf')


# In[80]:


tested.loc[tested['ORF'] == 'YOLO57W','ORF'] = 'YOL057W'
tested.loc[tested['ORF'] == 'YOLO62C','ORF'] = 'YOL062C'
tested.loc[tested['ORF'] == 'YJL206-A','ORF'] = 'YJL206C'
tested.loc[tested['ORF'] == 'YLR287-A','ORF'] = 'YLR287C-A'
tested.loc[tested['ORF'] == 'YBRF182C-A','ORF'] = 'YBR182C-A'


# In[81]:


# Make sure everything translated ok
t = looks_like_orf(tested['ORF'])
print(tested.loc[~t,])


# In[82]:


tested = tested.loc[t,:]


# In[84]:


tested = tested.drop_duplicates()


# In[113]:


missing = [orf for orf in data2.index.values if orf not in tested['ORF'].values]


# In[117]:


# In this case, the missing orf are DAMP strains which were included in the results but which we'll have to remove (only non-essential gene knockouts are kept)
data2 = data2.drop(index=missing)


# In[118]:


data2.shape


# # Prepare the final dataset

# In[119]:


dataset_ids = [16546,16590]


# In[120]:


datasets = datasets.reindex(index=dataset_ids)


# In[121]:


data = pd.DataFrame(index=tested['ORF'].values, columns=datasets['name'].values, data=0)


# In[122]:


data.loc[data2.index, :] = data2.values


# In[124]:


data.sum(axis=0)


# In[125]:


data = data.groupby(data.index).mean()


# In[126]:


# Create row index
data.index.name='orf'


# In[127]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[128]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[129]:


from IO.save_data_to_db2 import *


# In[130]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[131]:


save_data_to_db(data, paper_pmid)


# In[ ]:




