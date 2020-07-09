#!/usr/bin/env python
# coding: utf-8

# In[97]:


import numpy as np
import pandas as pd

from itertools import compress

import sys

from os.path import expanduser
sys.path.append(expanduser('~') + '/Lab/Utils/Python/')

from Conversions.translate import *
from Strings.is_a import *


# # Initial setup

# In[2]:


paper_pmid = 22384317
paper_name = 'fell_rosenwald_2011' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[24]:


original_data = pd.read_excel('raw_data/TableS1.xlsx', sheet_name='Sheet1', skiprows=1)


# In[25]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[26]:


original_data.head()


# In[27]:


original_data = original_data[['ORF','YPAD','HB']]


# In[28]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[29]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[30]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'], to='orf')


# In[31]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[70]:


phenotypes = sorted(pd.concat([original_data['YPAD'], original_data['HB']]).unique())


# In[71]:


phenotypes_dict = {k: 0 for k in phenotypes}


# In[72]:


phenotypes_dict['++++'] = 6
phenotypes_dict['+++'] = 5
phenotypes_dict['++'] = 4
phenotypes_dict['+'] = 3
phenotypes_dict['+/‐'] = 2
phenotypes_dict['‐'] = 1


# In[73]:


phenotypes_dict


# In[74]:


original_data['YPAD2'] = [phenotypes_dict[x.replace(' ','')] for x in original_data['YPAD'].values]


# In[75]:


original_data['HB2'] = [phenotypes_dict[x.replace(' ','')] for x in original_data['HB'].values]


# In[77]:


original_data.set_index('ORF', inplace=True)


# In[83]:


original_data['YPAD3'] = original_data['YPAD2'] - original_data.loc['WT','YPAD2']
original_data['HB3'] = original_data['HB2'] - original_data.loc['WT','HB2']


# In[84]:


original_data['HB3-YPAD3'] = original_data['HB3']-original_data['YPAD3']


# In[85]:


original_data.head()


# In[142]:


original_data.drop(index='WT', inplace=True)


# # Load & process tested strains

# In[124]:


tested = pd.read_csv('raw_data/Homo_diploids_041902.txt', sep='\t', header=1)


# In[125]:


tested = tested['ORF name'].unique()


# In[126]:


tested = tested.astype(str)


# In[127]:


tested = clean_orf(tested)


# In[128]:


tested = translate_sc(tested, to='orf')


# In[129]:


tested[tested=='YELOO1C'] = 'YEL001C'


# In[130]:


# Make sure everything translated ok
t = looks_like_orf(tested)


# In[131]:


list(compress(tested, ~np.array(t)))


# In[133]:


tested = [t for t in tested if t not in ['YMR41W','NAN']]


# In[135]:


# Test if all hits are present in tested
missing = [orf for orf in original_data.index.values if orf not in tested]
print(missing)


# In[136]:


tested = tested + ['YMR231W']


# # Prepare the final dataset

# In[145]:


dataset_ids = [16486]


# In[146]:


datasets = datasets.reindex(index=dataset_ids)


# In[147]:


data = pd.DataFrame(index=tested, columns=datasets['name'].values, data=0)


# In[148]:


data.loc[original_data.index, datasets['name'].values[0]] = original_data['HB3-YPAD3']


# In[149]:


data = data.groupby(data.index).mean()


# In[150]:


# Create row index
data.index.name='orf'


# In[151]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Print out

# In[153]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[154]:


from IO.save_data_to_db2 import *


# In[155]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[156]:


save_data_to_db(data, paper_pmid)


# In[ ]:




