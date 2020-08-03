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


paper_pmid = 17846143
paper_name = 'morton_coote_2007' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/Table1.xlsx', sheet_name='Sheet1')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[33]:


cols


# In[35]:


d_genes = []
mag_genes = []

cols = original_data.columns.values[2:]
for r in original_data.iterrows():
    for c in cols:
        s = str(r[1][c])
        s = s.replace('\xa0','')
        genes = s.split(',')
        if not isinstance(genes, list):
            genes = [genes]
        if c == cols[0]:
            d_genes = d_genes + genes
            mag_genes = mag_genes + genes
        elif c == cols[1]:
            d_genes = d_genes + genes
        elif c == cols[2]:
            mag_genes = mag_genes + genes


# In[39]:


d_genes = [s.strip() for s in d_genes if not s == 'nan']
mag_genes = [s.strip() for s in mag_genes if not s == 'nan']


# In[45]:


d_genes = clean_genename(d_genes)
mag_genes = clean_genename(mag_genes)


# In[46]:


d_orfs = translate_sc(d_genes, to='orf')
mag_orfs = translate_sc(mag_genes, to='orf')


# In[51]:


d_orfs = np.array(d_orfs)
mag_orfs = np.array(mag_orfs)


# In[56]:


d_orfs[d_orfs=='TMA29'] = 'YMR226C'


# In[57]:


t = looks_like_orf(d_orfs)
print(d_orfs[~np.array(t)])


# In[58]:


t = looks_like_orf(mag_orfs)
print(mag_orfs[~np.array(t)])


# In[60]:


all_orfs = np.unique(np.concatenate((d_orfs, mag_orfs)))


# In[63]:


data = pd.DataFrame(index=all_orfs, columns=['D','M'], data=np.zeros((len(all_orfs),2)))


# In[64]:


data.loc[d_orfs,'D'] = -1
data.loc[mag_orfs,'M'] = -1


# # Prepare the final dataset

# In[67]:


dataset_ids = [16536,16535]


# In[68]:


datasets = datasets.reindex(index=dataset_ids)


# In[69]:


data.columns = datasets['name'].values


# In[70]:


data = data.groupby(data.index).mean()


# In[71]:


# Create row index
data.index.name='orf'


# In[72]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[75]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[76]:


from IO.save_data_to_db2 import *


# In[77]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[79]:


save_data_to_db(data, paper_pmid)


# In[ ]:




