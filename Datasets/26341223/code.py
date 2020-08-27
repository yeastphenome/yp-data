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


paper_pmid = 26341223
paper_name = 'garcia_arroyo_2015' 


# In[232]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[233]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# ### Zymolyase

# In[5]:


path = 'raw_data/screening zymo/'
excel_files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]


# In[28]:


files_to_sheets = pd.read_excel('raw_data/zymo_files_to_sheets.xlsx', sheet_name='Sheet1', header=None)


# In[34]:


all_data = pd.DataFrame()
for f in excel_files:
    sheet_name = files_to_sheets.loc[files_to_sheets[0]==f,1].values[0]
    d = pd.read_excel(os.path.join(path, f), sheet_name=sheet_name)
    d.columns = [c.lower() for c in d.columns]
    all_data = pd.concat([all_data, d], axis=0)


# In[40]:


all_data.drop(columns=['unnamed: 3'], inplace=True)


# In[51]:


all_data['orf'] = all_data['orf'].astype(str)


# In[53]:


# Eliminate all white spaces & capitalize
all_data['orf'] = clean_orf(all_data['orf'])


# In[55]:


all_data = all_data.groupby('orf').mean()


# In[56]:


all_data = all_data.reset_index()


# In[57]:


# Translate to ORFs 
all_data['orf'] = translate_sc(all_data['orf'], to='orf')


# In[58]:


# Make sure everything translated ok
t = looks_like_orf(all_data['orf'])
print(all_data.loc[~t,])


# In[59]:


all_data.set_index('orf', inplace=True)


# In[60]:


all_data = all_data.div(all_data.loc['WT',:])


# In[62]:


all_data.drop(index='WT', inplace=True)


# In[127]:


all_data = all_data['24h'].to_frame()
all_data.columns = ['ZYM sens']


# In[190]:


all_data = all_data.groupby(all_data.index).mean()


# In[191]:


all_data.shape


# ### Congo Red and caspofungin

# In[63]:


all_data2 = pd.read_excel('raw_data/12864_2015_1879_MOESM2_ESM-2.xls', sheet_name='Hoja1', skiprows=1)


# In[66]:


all_data2['CR_num'] = all_data2['CR'].apply(lambda x: len(x) if isinstance(x, str) else x)


# In[68]:


all_data2['CAS_num'] = all_data2['CAS'].apply(lambda x: len(x) if isinstance(x, str) else x)


# In[69]:


all_data2['ZYM_num'] = all_data2['ZYM'].apply(lambda x: len(x) if isinstance(x, str) else x)


# In[71]:


all_data2 = all_data2[['ORF','CR_num','CAS_num','ZYM_num']]


# In[74]:


all_data2['ORF'] = all_data2['ORF'].astype(str)


# In[75]:


# Eliminate all white spaces & capitalize
all_data2['ORF'] = clean_orf(all_data2['ORF'])


# In[76]:


# Translate to ORFs 
all_data2['ORF'] = translate_sc(all_data2['ORF'], to='orf')


# In[77]:


# Make sure everything translated ok
t = looks_like_orf(all_data2['ORF'])
print(all_data2.loc[~t,])


# In[79]:


all_data2 = all_data2.loc[t,]


# In[82]:


all_data2 = all_data2.groupby('ORF').mean()


# In[85]:


all_data2 = -all_data2


# ### Caspofungin resistance

# In[155]:


all_data3 = pd.read_excel('raw_data/Table3.xlsx', sheet_name='Sheet1')


# In[156]:


all_data3['ORF'] = all_data3['ORF'].astype(str)


# In[157]:


# Eliminate all white spaces & capitalize
all_data3['ORF'] = clean_orf(all_data3['ORF'])


# In[158]:


# Translate to ORFs 
all_data3['ORF'] = translate_sc(all_data3['ORF'], to='orf')


# In[159]:


# Make sure everything translated ok
t = looks_like_orf(all_data3['ORF'])
print(all_data3.loc[~t,])


# In[160]:


all_data3 = all_data3[['ORF']].copy()


# In[162]:


all_data3['CAS res'] = 1


# In[163]:


all_data3 = all_data3.groupby('ORF').mean()


# # Prepare the final dataset

# In[234]:


data = all_data.join(all_data2, how='outer')


# In[235]:


data = data.join(all_data3, how='outer')


# In[236]:


data.head()


# In[237]:


# Dropping the discrete phenotypes for zymolyase because it is replaced by the quantitative data in ZYM sens
data.drop(columns=['ZYM_num'], inplace=True)


# In[238]:


dataset_ids = [16465,16466,16467,16470]


# In[239]:


datasets = datasets.reindex(index=dataset_ids)


# In[240]:


data.columns = datasets['name'].values


# In[241]:


data = data.groupby(data.index).mean()


# In[242]:


# Create row index
data.index.name='orf'


# In[243]:


print('Final data dimensions: %d x %d' % (data.shape))


# In[244]:


data.loc[:,col]


# In[245]:


# Set all NaN values to 0, effectively assuming that the quantitative ZYM sens dataset represents the tested universe for the other experiments as well
for c in [1,2,3]:
    col = data.columns.values[c]
    data.loc[data.loc[:,col].isnull(), col] = 0


# # Print out

# In[248]:


data.to_csv(paper_name + '.txt', sep='\t')


# # Save to DB

# In[249]:


from IO.save_data_to_db2 import *


# In[250]:


# Create column index
lst = [datasets.index.values, datasets['name'].values]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','dataset_name'])
data.columns = idx


# In[251]:


save_data_to_db(data, paper_pmid)


# In[ ]:




