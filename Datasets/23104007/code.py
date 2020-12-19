#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 23104007
paper_name = 'armakola_gitler_2012' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data 1

# In[5]:


original_data1 = pd.read_excel('raw_data/41588_2012_BFng2434_MOESM2_ESM.xlsx', sheet_name='Sys_Stnd_IDs_S_score_All.txt')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data1.shape))


# In[8]:


original_data1['orfs'] = original_data1['Systematic ID'].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data1['orfs'] = clean_orf(original_data1['orfs'])


# In[10]:


# Translate to ORFs 
original_data1['orfs'] = translate_sc(original_data1['orfs'], to='orf')


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data1['orfs'])
print(original_data1.loc[~t,])


# In[12]:


original_data1 = original_data1.loc[t,:]


# In[13]:


original_data1['data'] = original_data1['D S-score']


# In[14]:


original_data1.set_index('orfs', inplace=True)
original_data1.index.name='orf'


# In[15]:


original_data1 = original_data1.groupby(original_data1.index).mean()


# In[17]:


# Remove essentials (damp strains)
is_ess = is_essential(original_data1.index.values)


# In[20]:


original_data1 = original_data1.loc[~is_ess.values,:]


# In[21]:


original_data1.shape


# # Load & process data 2

# In[52]:


xl = pd.ExcelFile('raw_data/TDP-43 SSL complete round analysis rounds 1,2 8-16-08(norm version).xls')


# In[53]:


fields = ['ORF name','%D size-1 (rd 1)','%D size-2 (rd 1)','%D size-1 (rd 2)','%D size-2 (rd 2)']


# In[54]:


original_data2 = pd.DataFrame()
for sheet_name in xl.sheet_names:
    t = xl.parse(sheet_name)
    t = t.loc[:,fields]
    original_data2 = pd.concat((original_data2, t), axis=0)


# In[55]:


original_data2['orfs'] = original_data2['ORF name'].astype(str)


# In[56]:


# Eliminate all white spaces & capitalize
original_data2['orfs'] = clean_orf(original_data2['orfs'])


# In[57]:


original_data2 = original_data2.groupby('orfs').mean().reset_index()


# In[58]:


# Translate to ORFs 
original_data2['orfs'] = translate_sc(original_data2['orfs'], to='orf')


# In[60]:


original_data2.loc[original_data2['orfs']=='YLR287-A','orfs'] = 'YLR287C-A'


# In[61]:


# Make sure everything translated ok
t = looks_like_orf(original_data2['orfs'])
print(original_data2.loc[~t,])


# In[62]:


original_data2 = original_data2.loc[t,:]


# In[64]:


original_data2['data'] = original_data2[fields[1:]].mean(axis=1)


# In[65]:


original_data2.set_index('orfs', inplace=True)
original_data2.index.name='orf'


# In[66]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[70]:


original_data2.shape


# # Join the 2 datasets

# In[72]:


original_data = original_data1[['data']].join(original_data2[['data']], how='outer', lsuffix='_1', rsuffix='_2')


# In[73]:


original_data.head()


# # Prepare the final dataset

# In[74]:


data = original_data[['data_1','data_2']].copy()


# In[75]:


dataset_ids = [16494,16493]
datasets = datasets.reindex(index=dataset_ids)


# In[76]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[77]:


data.head()


# ## Subset to the genes currently in SGD

# In[78]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=original_data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[79]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[80]:


data.head()


# # Normalize

# In[81]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[82]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[83]:


data_norm[data.isnull()] = np.nan


# In[84]:


data_all = data.join(data_norm)


# In[85]:


data_all.head()


# # Print out

# In[86]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[87]:


from IO.save_data_to_db3 import *


# In[88]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




