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


paper_pmid = 23959798
paper_name = 'dong_rutherford_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/zmb999100152sd1.xlsx', sheet_name='Sheet1', skiprows=1)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[10]:


all_genes = []
for gene_list in original_data['Gene(s)'].values:
    all_genes = all_genes + gene_list.split(',')


# In[12]:


all_genes = list(set(all_genes))


# In[13]:


len(all_genes)


# In[15]:


# Eliminate all white spaces & capitalize
all_genes = clean_genename(all_genes)


# In[16]:


# Translate to ORFs 
all_orfs = translate_sc(all_genes, to='orf')


# In[18]:


# Make sure everything translated ok
t = looks_like_orf(all_orfs)
print(np.array(all_orfs)[~np.array(t)])


# # Load & process tested strains

# In[19]:


tested = pd.read_excel('raw_data/Gene list.xls', sheet_name='Sheet1', header=None)


# In[21]:


tested[1] = clean_orf(tested[1])


# In[22]:


tested[1] = translate_sc(tested[1], to='orf')


# In[23]:


# Make sure everything translated ok
t = looks_like_orf(tested[1])
print(tested.loc[~t])


# In[24]:


tested = tested[1].unique()


# In[27]:


missing = [orf for orf in all_orfs if orf not in tested]


# In[28]:


missing


# # Prepare the final dataset

# In[29]:


dataset_ids = [16559]


# In[30]:


datasets = datasets.reindex(index=dataset_ids)


# In[31]:


data = pd.DataFrame(index=tested, columns=datasets['name'].values, data=0)


# In[32]:


data.loc[all_orfs, datasets['name'].values[0]] = 1


# In[33]:


data = data.groupby(data.index).mean()


# In[34]:


# Create row index
data.index.name='orf'


# In[35]:


print('Final data dimensions: %d x %d' % (data.shape))


# In[36]:


data.sum(axis=0)


# # Print out

# In[37]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[38]:


from IO.save_data_to_db2 import *


# In[39]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[40]:


save_data_to_db(data, paper_pmid)


# In[ ]:




