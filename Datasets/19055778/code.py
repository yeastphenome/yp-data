#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 19055778
paper_name = 'alamgir_golshani_2008' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[63]:


original_data = pd.read_excel('raw_data/Plate_Analyzer_Alamgir_13mgmLParomomycin_08Dec05(1).xlsx', sheet_name='Sheet1')


# In[64]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[65]:


original_data.head()


# In[66]:


original_data.columns = ['plate','row','col','data']


# In[67]:


for c in ['plate','row','col']:
    original_data[c] = original_data[c].astype(int)


# In[68]:


print(original_data['plate'].unique())
print(original_data['row'].unique())
print(original_data['col'].unique())


# In[69]:


original_data = original_data.groupby(['plate','row','col']).mean()


# In[71]:


original_data.shape


# In[72]:


# Load platemap
files = ['Master plate list -- part 1.xlsx','Master plate list -- part2.xlsx']
mp_list = []
for f in files:
    mp = pd.read_excel('raw_data/' + f, sheet_name='Sheet2-Corrected')
#     print(mp.head())
    
    mp = mp[['Plate Number','Plate X','Plate Y','Systematic Name']].copy()
    mp.columns = ['plate','col','row','orf']
    
    mp_list.append(mp)


# In[73]:


mp = pd.concat(mp_list, axis=0, ignore_index=True)


# In[74]:


mp = mp.loc[mp.notnull().sum(axis=1)==4,:]


# In[75]:


for c in ['plate','row','col']:
    mp[c] = mp[c].astype(int)


# In[76]:


print(mp['plate'].unique())
print(mp['row'].unique())
print(mp['col'].unique())


# In[77]:


mp_count = mp.groupby(['plate','row','col']).size()


# In[79]:


(mp_count > 1).sum()


# In[80]:


mp.set_index(['plate','row','col'], inplace=True)


# In[81]:


original_data2 = original_data.join(mp, how='left')


# In[82]:


original_data2.head()


# In[83]:


original_data2.shape


# In[84]:


original_data2['orf'] = original_data2['orf'].astype(str)


# In[85]:


# Eliminate all white spaces & capitalize
original_data2['orf'] = clean_orf(original_data2['orf'])


# In[88]:


original_data2.loc[original_data2['orf'] == 'YPL072WA','orf'] = 'YPL072W-A'


# In[89]:


# Translate to ORFs 
original_data2['orf'] = translate_sc(original_data2['orf'].values, to='orf')


# In[90]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[91]:


original_data2.loc[~t & (np.abs(original_data2['data']) > 0)].shape


# In[92]:


original_data2 = original_data2.loc[t,:]


# In[93]:


original_data2['data'] = original_data2['data'].astype(float)


# In[94]:


original_data2.set_index('orf', inplace=True)


# In[95]:


original_data2 = original_data2[['data']].copy()


# In[96]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[97]:


original_data2.shape


# In[98]:


original_data2.head()


# # Prepare the final dataset

# In[105]:


data = original_data2.copy()


# In[106]:


dataset_ids = [142]
datasets = datasets.reindex(index=dataset_ids)


# In[107]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[108]:


data.head()


# ## Subset to the genes currently in SGD

# In[109]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[110]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[111]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[112]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[113]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[114]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[135]:


from IO.save_data_to_db3 import *


# In[136]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




