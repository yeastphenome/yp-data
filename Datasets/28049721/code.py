#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 28049721
paper_name = 'yofe_thoms_2017' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[57]:


original_data = pd.read_excel('raw_data/TablesS1-S7.xlsx', sheet_name='Table S1', skiprows=1)


# In[58]:


original_data.head()


# In[59]:


original_data['Deletion / DAmP'].unique()


# In[60]:


original_data = original_data.loc[original_data['Deletion / DAmP']=='Deletion',:].copy()


# In[61]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[62]:


original_data['orf'] = original_data['ORF'].astype(str)


# In[63]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[64]:


typo_fixes = {'YOLO57W':'YOL057W','YOLO62C':'YOL062C',
              'YKLO72W':'YKL072W','YJL206-A':'YJL206C-A','YLR287-A':'YLR287C-A'}


# In[65]:


for s in typo_fixes.keys():
    original_data.loc[original_data['orf']==s,'orf'] = typo_fixes[s]


# In[71]:


v[np.isnan(v)]


# In[72]:


# Translate to ORFs 
v = translate_sc(original_data['orf'].values, to='orf')
original_data['orf2'] = v


# In[73]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf2'])
print(original_data.loc[~t,])


# In[74]:


original_data.set_index('orf2', inplace=True)
original_data.index.name='orf'


# In[75]:


original_data = original_data[['Mean number of peroxisomes identified per cell [a.u]',
                              'Mean size of identified peroxisomes [a.u]']]


# In[76]:


original_data = original_data.groupby(original_data.index).mean()


# In[77]:


original_data.shape


# # Prepare the final dataset

# In[78]:


data = original_data.copy()


# In[79]:


dataset_ids = [11818, 11819]
datasets = datasets.reindex(index=dataset_ids)


# In[80]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[81]:


data.head()


# ## Subset to the genes currently in SGD

# In[82]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[83]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[84]:


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


# In[88]:


data_all.head()


# # Print out

# In[89]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[90]:


from IO.save_data_to_db3 import *


# In[91]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:



