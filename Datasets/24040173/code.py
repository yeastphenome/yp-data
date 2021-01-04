#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 24040173
paper_name = 'novo_gonzalez_2013' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[32]:


files = ['Only original foldchanges/Direct comparison/Phase I/HOP_t0vsYPD10 (dir) PhI.xlsx',
         'Only original foldchanges/Direct comparison/Phase I/HOP_t0vsMS10 (dir) PhI.xlsx',
         'Only original foldchanges/Direct comparison/Phase II/HOP_t0vsYPD10 (dir) PhII.xlsx',
         'Only original foldchanges/Direct comparison/Phase II/HOPt0vsHOP10 (dir) PhII.xlsx']

sheets = ['HOPt0vsYPD10_t0','HOPt0vsMS10_t0','HOPt0vsYPD10_t0','HOPt0vsHOP10_out0']


# In[33]:


original_data_list = []
for ixf, f in enumerate(files):
    original_data = pd.read_excel('raw_data/' + f, sheet_name=sheets[ixf])
    print('Original data dimensions: %d x %d' % (original_data.shape))
    
    original_data = original_data.loc[original_data['essential_gene']=='no',:]
    
    original_data['orf'] = original_data['strain'].apply(lambda x: x.split(':')[0])
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    # Make sure everything translated ok
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data = original_data.loc[t,:]
    
    original_data['data'] = original_data['Log2Ratio']
    original_data.set_index('orf', inplace=True)
    original_data = original_data[['data']].copy()
    original_data = original_data.groupby(original_data.index).mean()
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[34]:


original_data1 = pd.concat(original_data_list, axis=1)


# In[35]:


original_data1.columns = ['data1','data2','data3','data4']
original_data1.index.name='orf'


# In[36]:


original_data1['data1'] = original_data1['data1'] - original_data1['data2']
original_data1['data3'] = original_data1['data3'] - original_data1['data4']


# In[37]:


original_data1.drop(columns=['data2','data4'], inplace=True)


# In[38]:


original_data1.head()


# # Load & process data (het)

# In[23]:


files = ['Only original foldchanges/Direct comparison/Phase I/HIP_t0vsYPD20 (dir) PhI.xlsx',
         'Only original foldchanges/Direct comparison/Phase I/HIP_t0vsMS20 (dir) PhI.xlsx',
         'Only original foldchanges/Direct comparison/Phase II/HIP_t0vsYPD20 (dir) PhII.xlsx',
         'Only original foldchanges/Direct comparison/Phase II/HIPt0vsHIP10 (dir) PhII.xlsx']

sheets = ['HIPt0vsYPD20_t0','HIPt0vsMS20_tr0','HIPt0vsYPD20_t0','HIPt0vsHIP10_out0']


# In[24]:


original_data_list = []
for ixf, f in enumerate(files):
    original_data = pd.read_excel('raw_data/' + f, sheet_name=sheets[ixf])
    print('Original data dimensions: %d x %d' % (original_data.shape))
    
    original_data['orf'] = original_data['strain'].apply(lambda x: x.split(':')[0])
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    # Make sure everything translated ok
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data = original_data.loc[t,:]
    
    original_data['data'] = original_data['Log2Ratio']
    original_data.set_index('orf', inplace=True)
    original_data = original_data[['data']].copy()
    original_data = original_data.groupby(original_data.index).mean()
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[25]:


original_data2 = pd.concat(original_data_list, axis=1)


# In[26]:


original_data2.columns = ['data1','data2','data3','data4']
original_data2.index.name='orf'


# In[27]:


original_data2['data1'] = original_data2['data1'] - original_data2['data2']
original_data2['data3'] = original_data2['data3'] - original_data2['data4']


# In[28]:


original_data2.drop(columns=['data2','data4'], inplace=True)


# In[31]:


original_data2.head()


# # Merge all

# In[39]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[40]:


original_data.head()


# # Prepare the final dataset

# In[41]:


data = original_data.copy()


# In[42]:


dataset_ids = [197, 198, 15993, 15994]
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
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[46]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[47]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[48]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[49]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[50]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[51]:


from IO.save_data_to_db3 import *


# In[52]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




