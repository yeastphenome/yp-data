#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 15087496
paper_name = 'thorpe_dawes_2004' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[36]:


original_data = pd.read_excel('raw_data/05888Table2.xlsx', sheet_name='Total')


# In[37]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[38]:


original_data.head()


# In[39]:


original_data['orf'] = original_data['ORF'].astype(str)


# In[40]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[41]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[42]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[43]:


original_data.set_index('orf', inplace=True)


# In[44]:


original_data1 = original_data.loc[:,['CHP','Diamide','H2O2','LoaOOH','Menadione']].copy()
original_data2 = original_data.loc[:,['Unnamed: 4','Unnamed: 6','Unnamed: 8','Unnamed: 10','Unnamed: 12']].copy()


# In[45]:


data_switch = {'S': -1, 'R': 1}
for c in original_data1.columns:
    original_data1[c] = original_data1[c].apply(lambda x: data_switch[x] if x in data_switch.keys() else 0)


# In[46]:


original_data2 = original_data2.apply(pd.to_numeric, axis=1, errors='coerce')
original_data2 = original_data2 + 1
original_data2[original_data2.isnull()] = 1


# In[47]:


original_data = pd.DataFrame(index=original_data1.index, columns=original_data1.columns, data=np.multiply(original_data1.values, original_data2.values) )


# In[48]:


original_data.head()


# In[49]:


original_data = original_data.groupby(original_data.index).mean()


# In[50]:


original_data.shape


# # Prepare the final dataset

# In[51]:


data = original_data.copy()


# In[52]:


dataset_ids = [4959, 4957, 488, 4954, 4953]
datasets = datasets.reindex(index=dataset_ids)


# In[53]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[54]:


data.head()


# ## Subset to the genes currently in SGD

# In[55]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[56]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[57]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[58]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[59]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[60]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[61]:


from IO.save_data_to_db3 import *


# In[62]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:



