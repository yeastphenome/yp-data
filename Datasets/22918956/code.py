#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 22918956
paper_name = 'krumpe_schuldiner_2012' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[8]:


original_data = pd.read_csv('raw_data/hit_strains.txt', sep='\t')


# In[9]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[10]:


original_data.head()


# In[11]:


original_data['genes'] = original_data['Unnamed: 0'].astype(str)


# In[12]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[13]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['genes'], to='orf')


# In[14]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[15]:


original_data.set_index('orf', inplace=True)


# In[16]:


original_data = original_data.groupby(original_data.index).mean()


# In[17]:


original_data.shape


# # Load & process tested strains

# In[20]:


tested = pd.read_excel('raw_data/KO_DAmP_ORFs.xlsx', sheet_name='Sheet1', skiprows=1)


# In[21]:


tested.head()


# In[23]:


tested['orf'] = tested['ORF '].astype(str)


# In[24]:


tested['orf'] = clean_orf(tested['orf'])


# In[28]:


typo_fixes = {'YOLO57W':'YOL057W','YOLO62C':'YOL062C',
              'YBRF182C-A':'YBR182C-A','YLR287-A':'YLR287C-A','YJL206-A':'YJL206C-A'}
tested['orf'] = tested['orf'].apply(lambda x: typo_fixes[x] if x in typo_fixes.keys() else x)


# In[29]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[30]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[31]:


tested = tested.loc[t,:]


# In[32]:


tested_orfs = tested['orf'].unique()


# In[33]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[34]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[35]:


data = original_data.copy()


# In[36]:


dataset_ids = [15984, 16321]
datasets = datasets.reindex(index=dataset_ids)


# In[37]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[38]:


data.head()


# ## Subset to the genes currently in SGD

# In[39]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[40]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[41]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[42]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[43]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[44]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[45]:


from IO.save_data_to_db3 import *


# In[46]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




