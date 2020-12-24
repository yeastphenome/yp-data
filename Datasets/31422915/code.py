#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 31422915
paper_name = 'barbosa_siniossoglou_2019' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/mmc2.xlsx', sheet_name='Sheet1', skiprows=7)


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


d1 = 'ER+puncta'
d2 = 'Low signal'


# In[8]:


d1_ix = 0
d2_ix = original_data.loc[original_data.iloc[:,0]==d2].index.values[0]


# In[9]:


data1 = original_data.iloc[d1_ix:d2_ix,:].copy()
data2 = original_data.iloc[d2_ix+1:,:].copy()


# In[10]:


data1.columns = data1.iloc[0,:].copy()


# In[11]:


data2.columns = data2.iloc[0,:].copy()


# In[12]:


data1.head()


# In[13]:


c = 'Systematic Name'
data1[c] = data1[c].astype(str)
data2[c] = data2[c].astype(str)


# In[14]:


# Eliminate all white spaces & capitalize
data1[c] = clean_orf(data1[c])
data2[c] = clean_orf(data2[c])


# In[15]:


# Translate to ORFs 
data1[c] = translate_sc(data1[c], to='orf')
data2[c] = translate_sc(data2[c], to='orf')


# In[16]:


# Make sure everything translated ok
t = looks_like_orf(data1[c])
print(data1.loc[~t,])


# In[17]:


data1 = data1.loc[t,:]


# In[18]:


# Make sure everything translated ok
t = looks_like_orf(data2[c])
print(data2.loc[~t,])


# In[19]:


data2 = data2.loc[t,:]


# In[20]:


data1['data'] = 1
data2['data'] = 1


# In[21]:


data1.head()


# In[22]:


data1.set_index(c, inplace=True)
data2.set_index(c, inplace=True)


# In[23]:


data2.head()


# In[24]:


data = data1[['data']].join(data2[['data']], lsuffix='_1', rsuffix='_2', how='outer')


# In[25]:


data.loc[data['data_1'].isnull(),'data_1'] = 0
data.loc[data['data_2'].isnull(),'data_2'] = 0


# In[26]:


data.shape


# In[27]:


data.sum(axis=0)


# In[30]:


data.index.name='orf'


# In[31]:


data = data.groupby(data.index).mean()


# In[32]:


data2 = data.copy()


# In[33]:


data2.head()


# # Load & process tested strains

# In[34]:


tested = pd.read_excel('raw_data/KO_DAmP_ORFs.xlsx', sheet_name='Sheet1', skiprows=1)


# In[35]:


tested = tested.iloc[:,0].to_frame()


# In[36]:


tested.columns = ['ORF']


# In[37]:


tested['ORF'] = clean_orf(tested['ORF'])


# In[38]:


tested['ORF'] = translate_sc(tested['ORF'], to='orf')


# In[39]:


tested.loc[tested['ORF'] == 'YOLO57W','ORF'] = 'YOL057W'
tested.loc[tested['ORF'] == 'YOLO62C','ORF'] = 'YOL062C'
tested.loc[tested['ORF'] == 'YJL206-A','ORF'] = 'YJL206C'
tested.loc[tested['ORF'] == 'YLR287-A','ORF'] = 'YLR287C-A'
tested.loc[tested['ORF'] == 'YBRF182C-A','ORF'] = 'YBR182C-A'


# In[40]:


# Make sure everything translated ok
t = looks_like_orf(tested['ORF'])
print(tested.loc[~t,])


# In[41]:


tested = tested.loc[t,:]


# In[42]:


tested = tested.drop_duplicates()


# In[43]:


missing = [orf for orf in data2.index.values if orf not in tested['ORF'].values]


# In[44]:


missing


# In[46]:


# In this case, the missing orf are DAMP strains which were included in the results but which we'll have to remove (only non-essential gene knockouts are kept)
data2 = data2.drop(index=missing)


# In[48]:


data2.shape


# # Prepare the final dataset

# In[49]:


dataset_ids = [16546,16590]


# In[50]:


datasets = datasets.reindex(index=dataset_ids)


# In[51]:


data = data2.reindex(index=tested['ORF'].values, fill_value=0)


# In[52]:


data.sum(axis=0)


# In[53]:


data = data.groupby(data.index).mean()


# In[54]:


# Create row index
data.index.name='orf'


# In[56]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[57]:


print('Final data dimensions: %d x %d' % (data.shape))


# In[58]:


data.head()


# ## Subset to the genes currently in SGD

# In[59]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[60]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[61]:


data.head()


# # Normalize

# In[62]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[63]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[64]:


data_norm[data.isnull()] = np.nan


# In[65]:


data_all = data.join(data_norm)


# In[66]:


data_all.head()


# # Print out

# In[67]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[68]:


from IO.save_data_to_db3 import *


# In[69]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




