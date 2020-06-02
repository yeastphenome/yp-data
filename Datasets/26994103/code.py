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


paper_pmid = 26994103
paper_name = 'luo_jiang_2016' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[48]:


original_data1 = pd.read_excel('raw_data/Supplementary Table 2.xlsx', sheet_name='Sheet1', skiprows=2, names=['orf','gene','c1','c2'])
original_data2 = pd.read_excel('raw_data/Supplementary Table 3.xlsx', sheet_name='Sheet1', skiprows=2, names=['orf','gene','c1'])


# In[49]:


print('Original data dimensions: %d x %d' % (original_data1.shape))
print('Original data dimensions: %d x %d' % (original_data2.shape))


# In[50]:


original_data1['orf'] = original_data1['orf'].astype(str)
original_data2['orf'] = original_data2['orf'].astype(str)


# In[51]:


# Eliminate all white spaces & capitalize
original_data1['orf'] = clean_orf(original_data1['orf'])
original_data2['orf'] = clean_orf(original_data2['orf'])


# In[52]:


# Translate to ORFs 
original_data1['orf'] = translate_sc(original_data1['orf'], to='orf')
original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')


# In[53]:


# Make sure everything translated ok
t1 = looks_like_orf(original_data1['orf'])
t2 = looks_like_orf(original_data2['orf'])


# In[54]:


print(original_data1.loc[~t1,])


# In[55]:


print(original_data2.loc[~t2,])


# In[56]:


original_data1.set_index('orf', inplace=True)
original_data2.set_index('orf', inplace=True)


# In[57]:


original_data1['data'] = -1
original_data2['data'] = 1


# In[58]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[59]:


original_data = original_data[['data_1','data_2']].copy()


# In[61]:


original_data[original_data.isnull()] = 0


# In[62]:


original_data.head()


# In[65]:


# If the same strain is present more than once, average its values
data = original_data.groupby('orf').mean()


# In[66]:


print('Final data dimensions: %d x %d' % (data.shape))


# # Load the tested strains

# In[79]:


tested = pd.read_excel('raw_data/DELETION LIBRARY.xlsx', sheet_name='DELETION LIBRARY', skiprows=1)


# In[81]:


tested = tested['ORF name'].unique()


# In[89]:


tested = translate_sc(tested, to='orf')


# # Put together data and tested

# # Prepare the final dataset

# In[67]:


dataset_ids = [16448,16449]
datasets = datasets.reindex(index=dataset_ids)


# In[90]:


data = pd.DataFrame(index=tested, columns=datasets['name'].values, data=0)


# In[91]:


data.loc[original_data.index, datasets['name'].values[0]] = original_data['data_1']


# In[92]:


data.loc[original_data.index, datasets['name'].values[1]] = original_data['data_2']


# In[97]:


# Create row index
data.index.name='orf'


# # Print out

# In[98]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[99]:


from IO.save_data_to_db2 import *


# In[104]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[106]:


save_data_to_db(data, paper_pmid)


# In[107]:


data.head()


# In[109]:


datasets = data.columns.values
orfs = data.index.values


# In[111]:


pmid = paper_pmid


# In[113]:


sql_del = 'DELETE FROM datasets_data WHERE dataset_id in '           '(SELECT datasets_dataset.id FROM datasets_dataset '           'INNER JOIN papers_paper ON datasets_dataset.paper_id = papers_paper.id '           'WHERE pmid = ' + str(pmid) + ');\n'

sql_in = 'BEGIN;\n'
sql_template = 'INSERT INTO datasets_data (dataset_id, orf, value) VALUES (%d, \'%s\', %f);\n'

for dataset in datasets:
    for orf in orfs:
        if not np.isnan(data.loc[orf, dataset]):
            sql_in += (sql_template % (dataset, orf, data.loc[orf, dataset]))


# In[ ]:





# In[116]:


data.loc[orf, dataset]


# In[ ]:




