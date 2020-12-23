#!/usr/bin/env python
# coding: utf-8

# In[3]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[4]:


paper_pmid = 31885205
paper_name = 'galardini_beltrao_2019' 


# In[5]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[6]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[12]:


original_data = pd.read_csv('raw_data/Table_EV2_ko_scores.txt', sep='\t', skiprows=2)


# In[15]:


original_data = original_data.loc[original_data['strain']=='S288C',:]


# In[16]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[21]:


original_data2 = original_data.groupby(['condition','gene'])['score'].mean().reset_index()


# In[23]:


original_data3 = pd.pivot_table(original_data2, index='gene', columns='condition', values='score')


# In[24]:


original_data3.shape


# In[33]:


original_data3['orfs'] = original_data3.index.values.astype(str)


# In[34]:


# Eliminate all white spaces & capitalize
original_data3['orfs'] = clean_orf(original_data3['orfs'])


# In[35]:


# Fix typos
original_data3.loc[original_data3['orfs']=='YLR287-A','orfs'] = 'YLR287C-A'


# In[36]:


# Translate to ORFs 
original_data3['orfs'] = translate_sc(original_data3['orfs'], to='orf')


# In[37]:


# Make sure everything translated ok
t = looks_like_orf(original_data3['orfs'])
print(original_data3.loc[~t,])
# Output: found 22 inexisting ORFs and WT


# In[38]:


original_data3 = original_data3.loc[t,:]


# In[39]:


original_data3.set_index('orfs', inplace=True)
original_data3.index.name='orf'


# In[40]:


original_data3 = original_data3.groupby(original_data3.index).mean()


# In[41]:


original_data3.shape


# # Load dataset ids

# In[62]:


dataset_id_map = pd.read_excel('raw_data/msb198831-sup-0002-tableev1.xlsx', encoding='utf-8',
                               sheet_name='Sheet1', header=None)


# In[63]:


dataset_id_map[0] = dataset_id_map[0].apply(lambda x: x.replace('Â',''))


# In[65]:


dataset_id_map.set_index(0, inplace=True)


# In[66]:


dataset_id_map = dataset_id_map.reindex(index=original_data3.columns.values)


# In[68]:


dataset_ids = list(dataset_id_map[1])


# # Prepare the final dataset

# In[70]:


data = original_data3.copy()


# In[71]:


datasets = datasets.reindex(index=dataset_ids)


# In[72]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[73]:


data.head()


# ## Subset to the genes currently in SGD

# In[74]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[75]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[76]:


data.head()


# # Normalize

# In[77]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[78]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[79]:


data_norm[data.isnull()] = np.nan


# In[80]:


data_all = data.join(data_norm)


# In[81]:


data_all.head()


# # Print out

# In[82]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[83]:


from IO.save_data_to_db3 import *


# In[84]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




