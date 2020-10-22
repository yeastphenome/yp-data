#!/usr/bin/env python
# coding: utf-8

# In[45]:


import numpy as np
import pandas as pd

import sys

from os.path import expanduser
sys.path.append(expanduser('~') + '/Lab/Utils/Python/')

from Conversions.translate import *
from Strings.is_a import *


# # Initial setup

# In[46]:


paper_pmid = 17630978
paper_name = 'pagani_arino_2007' 


# In[47]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[48]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[49]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, names=['genes','data'], sep='\t')


# In[50]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[51]:


original_data['genes'] = original_data['genes'].astype(str)


# In[52]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[53]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[54]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[55]:


original_data.set_index('orfs', inplace=True)


# # Load & process tested strains

# In[78]:


tested = pd.read_excel('raw_data/EUROFAN haploid collection .xlsx', sheet_name='GENERAL 07-02-11')


# In[79]:


tested['orf'] = tested['Systematic Name '].astype(str)
tested['orf'] = clean_orf(tested['orf'])


# In[80]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[81]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[82]:


tested.loc[~t,'orf'].unique()


# In[83]:


tested = tested.loc[t,:]


# In[84]:


tested = tested['orf'].unique()


# In[85]:


tested.shape


# In[86]:


missing = [orf for orf in original_data.index.values if orf not in tested]
missing


# In[88]:


# Decided to add the 1 missing strain
tested = np.append(tested, 'YBR011C')


# # Prepare the final dataset

# In[89]:


dataset_ids = [16619]


# In[90]:


datasets = datasets.reindex(index=dataset_ids)


# In[91]:


data = pd.DataFrame(index=tested, columns=datasets['name'].values, data=0)


# In[92]:


data.loc[original_data.index, datasets['name'].values[0]] = original_data['data']


# In[93]:


data = data.groupby(data.index).mean()


# In[94]:


# Create row index
data.index.name='orf'


# In[95]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[99]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[100]:


from IO.save_data_to_db2 import *


# In[101]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[102]:


save_data_to_db(data, paper_pmid)


# In[ ]:




