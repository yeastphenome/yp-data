#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 28152593
paper_name = 'lebesgue_lemeer_2017' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[8]:


original_data = pd.read_excel('raw_data/pr6b00691_si_001.xlsx', sheet_name='S5-Deletion-Overexpression ', skiprows=1)


# In[9]:


original_data.head()


# In[10]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[11]:


original_data['orf'] = original_data['ORF (SGD)'].astype(str)


# In[12]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[15]:


typo_fixes = {'YBRF182C-A':'YBR182C-A','YJL206-A':'YJL206C-A','YKLO72W':'YKL072W',
              'YLR287-A':'YLR287C-A','YOLO57W':'YOL057W','YOLO62C':'YOL062C'}


# In[16]:


for s in typo_fixes.keys():
    original_data.loc[original_data['orf']==s,'orf'] = typo_fixes[s]


# In[17]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[18]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[19]:


original_data['data'] = original_data['Size ratio after deletion'].astype(float)


# In[20]:


original_data.set_index('orf', inplace=True)


# In[21]:


original_data = original_data[['data']]


# In[22]:


original_data = original_data.groupby(original_data.index).mean()


# In[23]:


original_data.shape


# In[24]:


dataset_ids = [11864]


# # Prepare the final dataset

# In[25]:


data = original_data.copy()


# In[26]:


datasets = datasets.reindex(index=dataset_ids)


# In[27]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[28]:


data.head()


# ## Subset to the genes currently in SGD

# In[29]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[30]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[31]:


data.head()


# # Normalize

# In[32]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[33]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[34]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[35]:


data_all.head()


# # Print out

# In[36]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[37]:


from IO.save_data_to_db3 import *


# In[38]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:



