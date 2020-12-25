#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 28818866
paper_name = 'segura_wang_korbel_2017' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[23]:


original_data_all = pd.read_excel('raw_data/TableS5.xlsx', sheet_name='Tabelle1', skiprows=2)


# In[24]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[25]:


drugs = ['Campt','Doxo','HU','MMS']


# In[26]:


original_data_list = []
for di, d in enumerate(drugs):
    start = original_data_all.loc[original_data_all['Strain'] == drugs[di]].index.values[0]
    if di < 3:
        stop = original_data_all.loc[original_data_all['Strain'] == drugs[di+1]].index.values[0]
    else:
        stop = -1
    
    original_data = original_data_all.iloc[start:stop].copy()
    original_data['orf'] = original_data['Strain'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data = original_data.loc[t,:]
    
    original_data.set_index('orf', inplace=True)

    original_data['data'] = (original_data['Fold enrich. in Drug'] - original_data['Fold enrich. in YPAD'])/original_data['Fold enrich. in YPAD']
    original_data = original_data[['data']].copy()
    original_data = original_data.groupby(original_data.index).mean()
    
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[27]:


original_data1, original_data2, original_data3, original_data4 = original_data_list


# In[28]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')
original_data = original_data.join(original_data3, how='outer', lsuffix='', rsuffix='_3')
original_data = original_data.join(original_data4, how='outer', lsuffix='', rsuffix='_4')


# In[36]:


original_data[original_data.isnull()] = 0


# In[37]:


original_data.shape


# In[38]:


original_data.head()


# In[39]:


dataset_ids = [16150, 16149, 16147, 16148]


# # Prepare the final dataset

# In[40]:


data = original_data.copy()


# In[41]:


datasets = datasets.reindex(index=dataset_ids)


# In[42]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[43]:


data.head()


# ## Subset to the genes currently in SGD

# In[44]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[45]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[46]:


data.head()


# # Normalize

# In[47]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[48]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[49]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[50]:


data_all.head()


# # Print out

# In[51]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[52]:


from IO.save_data_to_db3 import *


# In[53]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




