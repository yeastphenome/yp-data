#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 28598353
paper_name = 'zatorska_strahl_2017' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[27]:


original_data = pd.read_excel('raw_data/ijms-18-01226-s001.xlsx', sheet_name='Table S1', skiprows=2)


# In[28]:


original_data.head()


# In[29]:


original_data = original_data.loc[original_data['Library']=='deletion',:]


# In[30]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[31]:


original_data['orf'] = original_data['Systematic name'].astype(str)


# In[32]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[33]:


typo_fix = {'YLR287-A':'YLR287C-A','YBRF182C-A':'YBR182C-A','YOLO57W':'YOL057W',
            'YOLO62C':'YOL062C','YBRF182C-A':'YBR182C-A','YKLO72W':'YKL072W',
           'YJL206-A':'YJL206C-A'}


# In[34]:


for s in typo_fix.keys():
    original_data.loc[original_data['orf']==s,'orf'] = typo_fix[s]


# In[35]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[36]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[37]:


original_data = original_data.loc[t,:]


# In[38]:


original_data.set_index('orf', inplace=True)


# In[40]:


original_data['data'] = original_data['Growth [%] +inhibitor/-inhibitor']


# In[41]:


original_data = original_data[['data']].copy()


# In[42]:


original_data = original_data.groupby(original_data.index).mean()


# In[43]:


original_data.shape


# # Prepare the final dataset

# In[44]:


data = original_data[['data']].copy()


# In[45]:


dataset_ids = [11862]
datasets = datasets.reindex(index=dataset_ids)


# In[46]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[47]:


data.head()


# ## Subset to the genes currently in SGD

# In[48]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[49]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[50]:


data.head()


# # Normalize

# In[51]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[52]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[53]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[54]:


data_all.head()


# # Print out

# In[55]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[56]:


from IO.save_data_to_db3 import *


# In[57]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:



