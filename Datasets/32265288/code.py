#!/usr/bin/env python
# coding: utf-8

# In[1]:


import numpy as np
import pandas as pd

import sys
import regex

from os.path import expanduser
sys.path.append(expanduser('~') + '/Lab/Utils/Python/')

from Conversions.translate import *
from Strings.is_a import *


# # Initial setup

# In[2]:


paper_pmid = 32265288
paper_name = 'novarina_chang_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data1 = pd.read_excel('raw_data/Supp_Tables_v4.xlsx', sheet_name='Table S2', skiprows=2)
original_data2 = pd.read_excel('raw_data/Supp_Tables_v4.xlsx', sheet_name='Table S4', skiprows=2)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))


# In[7]:


print('Original data dimensions: %d x %d' % (original_data2.shape))


# ### Dataset 1

# In[8]:


gene_col1 = 'Positives from patch assay'


# In[9]:


original_data1['genes'] = original_data1[gene_col1].astype(str)


# In[10]:


# Eliminate all white spaces & capitalize
original_data1['genes'] = clean_genename(original_data1['genes'])


# In[11]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['genes'], to='orf')


# In[12]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[13]:


original_data1 = original_data1.loc[t,:]


# In[14]:


original_data1['data'] = 1


# In[15]:


original_data1.set_index('orfs', inplace=True)


# In[49]:


original_data1 = original_data1.groupby(original_data1.index).mean()


# ### Dataset 2

# In[16]:


orf_col = 'ORF name'


# In[17]:


original_data2['orfs'] = original_data2[orf_col].astype(str)


# In[18]:


# Eliminate all white spaces & capitalize
original_data2['orfs'] = clean_orf(original_data2['orfs'])


# In[19]:


# Translate to ORFs 
original_data2['orfs'] = translate_sc(original_data2['orfs'], to='orf')


# In[20]:


original_data2.loc[original_data2['orfs']=='YLR287-A','orfs'] = 'YLR287C-A'


# In[21]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[22]:


wt_phenotype = 56 # from main text
original_data2['data'] = (original_data2['Percent recombinants'] / wt_phenotype) - 1


# In[23]:


original_data2.set_index('orfs', inplace=True)


# In[50]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# ### Join the 2 datasets

# In[51]:


original_data = original_data1[['data']].join(original_data2[['data']], how='outer', lsuffix='_1', rsuffix='_2')


# In[52]:


original_data.shape


# In[53]:


# Orfs from dataset 1 that were not present in dataset 2
unmatched_orfs = original_data.loc[original_data['data_2'].isnull()].index.values


# In[54]:


# Check for partial matches of hits from dataset1 in the list from dataset2 (possible typos in dataset1)
orfs_ref = pd.Series(original_data2.index.values)

for orf in unmatched_orfs:
    s = '(' + orf + '){e<=1}'
    
    m = orfs_ref.apply(lambda x: len(regex.findall(s, x)))
    nm = m.sum()
    
    print('%s\t%d' % (orf, nm))


# In[55]:


# No obvious typo fixes. Decided to leave dataset 2 values for unmatched orfs at NaN. But use dataset 2 list as the tested universe for dataset 1 (best approximation).
original_data['data_1'].loc[original_data['data_1'].isnull()] = 0


# In[56]:


original_data.notnull().sum(axis=0)


# # Prepare the final dataset

# In[57]:


dataset_ids = [16617, 16618]


# In[58]:


datasets = datasets.reindex(index=dataset_ids)


# In[59]:


data = original_data[['data_1','data_2']].copy()


# In[60]:


data.columns = datasets['name'].values


# In[61]:


data = data.groupby(data.index).mean()


# In[62]:


# Create row index
data.index.name='orf'


# In[63]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[64]:


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




