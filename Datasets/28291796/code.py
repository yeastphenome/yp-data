#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 28291796
paper_name = 'huseinovic_vos_2017' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[16]:


original_data1 = pd.read_excel('raw_data/pone.0173573.s003.xls', sheet_name='Integrated intensity', skiprows=5)


# In[17]:


print('Original data dimensions: %d x %d' % (original_data1.shape))


# In[18]:


original_data1.head()


# In[20]:


original_data1a = original_data1[['ORF','Integrated intensity','Integrated intensity.1']]
original_data1b = original_data1[['ORF.1','Integrated intensity.2','Integrated intensity.3']]
original_data1c = original_data1[['ORF.2','Integrated intensity.4','Integrated intensity.5']]
original_data1d = original_data1[['ORF.3','Integrated intensity.6','Integrated intensity.7']]

cols = ['orf','data1','data2']
original_data1a.columns = cols
original_data1b.columns = cols
original_data1c.columns = cols
original_data1d.columns = cols


# In[37]:


original_data1 = pd.concat([original_data1a, original_data1b, original_data1c, original_data1d], axis=0, ignore_index=True)


# In[38]:


original_data1.head()


# In[39]:


original_data1['orf'] = original_data1['orf'].astype(str)


# In[40]:


# Eliminate all white spaces & capitalize
original_data1['orf'] = clean_orf(original_data1['orf'])


# In[41]:


# Translate to ORFs 
original_data1['orf'] = translate_sc(original_data1['orf'], to='orf')


# In[42]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orf'])
print(original_data1.loc[~t,])


# In[43]:


original_data1 = original_data1.loc[t,]


# In[44]:


for d in ['data1','data2']:
    original_data1[d] = pd.to_numeric(original_data1[d], errors='coerce')


# In[45]:


original_data1['data'] = original_data1[['data1','data2']].mean(axis=1)


# In[46]:


original_data1.set_index('orf', inplace=True)


# In[47]:


original_data1 = original_data1[['data']].copy()


# In[48]:


original_data1 = original_data1.groupby(original_data1.index).mean()


# In[49]:


original_data1.shape


# In[50]:


original_data1.head()


# # Load & process the data (2)

# In[55]:


original_data2 = pd.read_excel('raw_data/Krogan screen all data row sizes 27sep16.xlsx', sheet_name='30C day2 ', skiprows=3)
original_data3 = pd.read_excel('raw_data/Krogan screen all data row sizes 27sep16.xlsx', sheet_name='37C day2', skiprows=3)


# In[58]:


print('Original data dimensions: %d x %d' % (original_data2.shape))
print('Original data dimensions: %d x %d' % (original_data3.shape))


# In[59]:


original_data2.head()


# In[60]:


original_data2['orf'] = original_data2['ORF'].astype(str)
original_data3['orf'] = original_data3['ORF'].astype(str)


# In[61]:


# Eliminate all white spaces & capitalize
original_data2['orf'] = clean_orf(original_data2['orf'])
original_data3['orf'] = clean_orf(original_data3['orf'])


# In[62]:


# Translate to ORFs 
original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')
original_data3['orf'] = translate_sc(original_data3['orf'], to='orf')


# In[63]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orf'])
print(original_data2.loc[~t,])


# In[64]:


# Make sure everything translated ok
t = looks_like_orf(original_data3['orf'])
print(original_data3.loc[~t,])


# In[65]:


original_data2 = original_data2.loc[t,]
original_data3 = original_data3.loc[t,]


# In[68]:


original_data3.head()


# In[67]:


for d in ['Raw colony size','Raw colony size.1','Raw colony size.2','Raw colony size.3','Raw colony size.4','Raw colony size.5','Raw colony size.6']:
    original_data2[d] = pd.to_numeric(original_data2[d], errors='coerce')


# In[69]:


for d in ['Raw colony size','Raw colony size.1','Raw colony size.2','Raw colony size.3','Raw colony size.4','Raw colony size.5','Raw colony size.6','Raw colony size.7','Raw colony size.8','Raw colony size.9']:
    original_data3[d] = pd.to_numeric(original_data3[d], errors='coerce')


# In[70]:


original_data3['0'] = original_data3[['Raw colony size','Raw colony size.1']].mean(axis=1)
original_data3['50'] = original_data3[['Raw colony size.2','Raw colony size.3']].mean(axis=1)
original_data3['60'] = original_data3[['Raw colony size.4','Raw colony size.5']].mean(axis=1)
original_data3['70'] = original_data3[['Raw colony size.6','Raw colony size.7']].mean(axis=1)
original_data3['80'] = original_data3[['Raw colony size.8','Raw colony size.9']].mean(axis=1)


# In[71]:


original_data2.set_index('orf', inplace=True)
original_data3.set_index('orf', inplace=True)


# In[74]:


original_data2 = original_data2[['Raw colony size','Raw colony size.1','Raw colony size.2','Raw colony size.3','Raw colony size.4','Raw colony size.5','Raw colony size.6']].copy()
original_data2.columns = ['0','50','60','70','80','90','100']


# In[75]:


original_data3 = original_data3[['0','50','60','70','80']].copy()


# In[76]:


original_data4 = original_data2.join(original_data3, how='outer', lsuffix='_30', rsuffix='_37')


# In[78]:


original_data4.shape


# In[79]:


original_data4 = original_data4.groupby(original_data4.index).mean()


# In[80]:


original_data4.shape


# In[81]:


original_data4.head()


# In[82]:


# Normalize to dose 0
for c in ['50_30','60_30','70_30','80_30','90','100']:
    original_data4[c] = original_data4[c] / original_data4['0_30']
    
for c in ['50_37','60_37','70_37','80_37']:
    original_data4[c] = original_data4[c] / original_data4['0_37']


# In[84]:


original_data4.drop(columns=['0_30','0_37'], inplace=True)


# In[85]:


original_data4.head()


# In[86]:


# Merge all datasets
original_data = original_data1.join(original_data4, how='outer')


# In[88]:


original_data.shape


# In[90]:


original_data.head()


# # Prepare the final dataset

# In[89]:


data = original_data.copy()


# In[91]:


dataset_ids = [22075, 22088, 22089, 22090, 22091, 22092, 22087, 22076, 22082, 22083, 22084]
datasets = datasets.reindex(index=dataset_ids)


# In[92]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[93]:


data.head()


# ## Subset to the genes currently in SGD

# In[94]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[95]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[96]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[97]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[98]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[99]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[100]:


from IO.save_data_to_db3 import *


# In[101]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




