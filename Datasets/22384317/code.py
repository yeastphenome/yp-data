#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 22384317
paper_name = 'fell_rosenwald_2011' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[45]:


original_data = pd.read_excel('raw_data/TableS1.xlsx', sheet_name='Sheet1', skiprows=1)


# In[46]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[47]:


original_data.head()


# In[48]:


original_data = original_data[['ORF','YPAD','HB']]


# In[49]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[50]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[51]:


# Translate to ORFs 
original_data['ORF'] = translate_sc(original_data['ORF'], to='orf')


# In[52]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])
print(original_data.loc[~t,])


# In[53]:


phenotypes = sorted(pd.concat([original_data['YPAD'], original_data['HB']]).unique())


# In[54]:


phenotypes_dict = {k: 0 for k in phenotypes}


# In[55]:


phenotypes_dict['++++'] = 6
phenotypes_dict['+++'] = 5
phenotypes_dict['++'] = 4
phenotypes_dict['+'] = 3
phenotypes_dict['+/‐'] = 2
phenotypes_dict['‐'] = 1


# In[56]:


phenotypes_dict


# In[57]:


original_data['YPAD2'] = [phenotypes_dict[x.replace(' ','')] for x in original_data['YPAD'].values]


# In[58]:


original_data['HB2'] = [phenotypes_dict[x.replace(' ','')] for x in original_data['HB'].values]


# In[59]:


original_data.set_index('ORF', inplace=True)
original_data.index.name='orf'


# In[60]:


original_data['YPAD3'] = original_data['YPAD2'] - original_data.loc['WT','YPAD2']
original_data['HB3'] = original_data['HB2'] - original_data.loc['WT','HB2']


# In[61]:


original_data['HB3-YPAD3'] = original_data['HB3']-original_data['YPAD3']


# In[62]:


original_data.head()


# In[63]:


original_data.drop(index='WT', inplace=True)


# In[64]:


original_data['data'] = original_data['HB3-YPAD3']
original_data = original_data[['data']].copy()


# In[65]:


original_data = original_data.groupby(original_data.index).mean()


# In[66]:


original_data.shape


# # Load & process tested strains

# In[67]:


tested = pd.read_csv('raw_data/Homo_diploids_041902.txt', sep='\t', header=1)


# In[68]:


tested['orf'] = tested['ORF name'].astype(str)


# In[69]:


tested['orf'] = clean_orf(tested['orf'])


# In[70]:


tested.loc[tested['orf']=='YELOO1C','orf'] = 'YEL001C'


# In[71]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[72]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[73]:


tested = tested.loc[t,:]


# In[74]:


tested_orfs = tested['orf'].unique()


# In[75]:


# Test if all hits are present in tested
missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
print(missing)


# In[76]:


tested_orfs = list(tested_orfs) + ['YMR231W']


# In[77]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# In[78]:


original_data.head()


# # Prepare the final dataset

# In[79]:


data = original_data.copy()


# In[80]:


dataset_ids = [16486]
datasets = datasets.reindex(index=dataset_ids)


# In[81]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[82]:


data.head()


# ## Subset to the genes currently in SGD

# In[83]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[84]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[85]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[86]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[87]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[88]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[89]:


from IO.save_data_to_db3 import *


# In[90]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




