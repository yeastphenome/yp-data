#!/usr/bin/env python
# coding: utf-8

# In[23]:


import numpy as np
import pandas as pd

import sys
import re

from os.path import expanduser
sys.path.append(expanduser('~') + '/Lab/Utils/Python/')

from Conversions.translate import *
from Strings.is_a import *


# # Initial setup

# In[4]:


paper_pmid = 18755837
paper_name = 'huang_bystrom_2008' 


# In[5]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[6]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[13]:


original_data = pd.read_csv('raw_data/Table1-2.txt', header=None, names=['genes','score'], sep='\t')


# In[14]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[15]:


original_data['genes'] = original_data['genes'].astype(str)


# In[16]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[17]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes'], to='orf')


# In[18]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[19]:


original_data = original_data.loc[t,:]


# In[20]:


original_data.set_index('orfs', inplace=True)


# # Load & process tested strains

# In[60]:


orf_pattern = 'Y[A-P][RL][0-9]{3}[CW](-[A-H])*'
p = re.compile(orf_pattern)


# In[64]:


file_path = 'raw_data/Deletion collection homo dipl.txt'

tested_orfs = []
with open(file_path,'r') as fp:
    line = fp.readline()
    while line:
        res = p.search(line)
        if res:
            tested_orfs.append(res.group(0))
        line = fp.readline()
        


# In[67]:


tested_orfs = np.unique(np.array(tested_orfs))


# In[69]:


tested_orfs = clean_orf(tested_orfs)


# In[70]:


tested_orfs2 = translate_sc(tested_orfs, to='orf')


# In[71]:


# Make sure everything translated ok
t = looks_like_orf(tested_orfs2)
print(tested_orfs2[~np.array(t)])


# In[79]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs2]


# In[81]:


missing


# # Prepare the final dataset

# In[82]:


dataset_ids = [16436]


# In[83]:


datasets = datasets.reindex(index=dataset_ids)


# In[84]:


data = pd.DataFrame(index=tested_orfs2, columns=datasets['name'].values, data=0)


# In[85]:


data.loc[original_data.index, datasets['name'].values[0]] = original_data['score']


# In[86]:


data = data.groupby(data.index).mean()


# In[87]:


# Create row index
data.index.name='orf'


# In[88]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[89]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[90]:


from IO.save_data_to_db2 import *


# In[91]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[92]:


save_data_to_db(data, paper_pmid)


# In[ ]:




