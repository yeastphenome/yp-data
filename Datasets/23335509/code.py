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
from Math.normalize import z_transform_mode


# # Initial setup

# In[2]:


paper_pmid = 23335509
paper_name = 'cuesta_marban_mollinedo_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# In[5]:


path_to_genes = '../../Private-Utils/datasets_gene2.txt'
path_to_consensus_tested = '../../Private-Utils/yp_2020-09-01_orfs.txt'


# # Load & process the data

# In[26]:


original_data = pd.read_excel('raw_data/jbc.M112.425769-1.xlsx', sheet_name='Sheet1')


# In[27]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[28]:


original_data['genes'] = original_data['Gene'].astype(str)


# In[29]:


original_data['genes2'] = original_data['genes']

ix = original_data['genes2'].str.contains('/')
for i, r in original_data.loc[ix,].iterrows():
    t = r['genes2'].split('/')
    original_data.at[i,'genes2'] = t[0]


# In[30]:


# Eliminate all white spaces & capitalize
original_data['genes2'] = clean_genename(original_data['genes2'])


# In[31]:


# Translate to ORFs 
original_data['orfs'] = translate_sc(original_data['genes2'], to='orf')


# In[32]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orfs'])
print(original_data.loc[~t,])


# In[33]:


original_data['data'] = original_data['Resistance'].apply(lambda x: len(x))


# In[38]:


original_data.set_index('orfs', inplace=True)
original_data.index.name='orf'


# In[39]:


original_data = original_data.groupby(original_data.index).mean()


# In[40]:


original_data.shape


# # Prepare the final dataset

# In[41]:


data = original_data[['data']].copy()


# In[42]:


dataset_ids = [16563]
datasets = datasets.reindex(index=dataset_ids)


# In[43]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[44]:


data.head()


# ## Subset to the genes currently in SGD

# In[45]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=original_data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[46]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[47]:


data.head()


# # Normalize
# 

# In[55]:


def normalize_phenotypic_scores(df, has_tested=False):
    
    if not has_tested:
        
        genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
        genes = genes.reset_index().set_index('systematic_name', drop=False)
        
        yp_orfs = pd.read_csv(path_to_consensus_tested, header=None)
        yp_orfs = yp_orfs[1].values
        consensus_tested = [tuple(genes.loc[orf,['id','systematic_name']]) for orf in yp_orfs]
        
        df_index = [tuple(x) for x in df.index]
        
        consensus_tested = list(set(consensus_tested + df_index))
        
        df = df.reindex(index=consensus_tested, fill_value=0)
        
    df_norm = z_transform_mode(df)
    
    return df_norm


# In[56]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[57]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[58]:


data_norm[data.isnull()] = np.nan


# In[59]:


data_all = data.join(data_norm)


# In[60]:


data_all.head()


# # Print out

# In[62]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[63]:


from IO.save_data_to_db3 import *


# In[64]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




