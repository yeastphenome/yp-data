#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 24968945
paper_name = 'gaupel_tenniswood_2004' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


sheets = ['Top sensitive strains','Top resistant strains']


# In[7]:


original_data_list = []
for s in sheets:
    original_data = pd.read_excel('raw_data/1471-2164-15-528-s1.xlsx', sheet_name=s)
    print('Original data dimensions: %d x %d' % (original_data.shape))
    original_data['orf'] = original_data['Unique Name'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    
    original_data.loc[original_data['orf']=='YJL206-A','orf'] = 'YJL206C-A'
    
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    
    # Make sure everything translated ok
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    
    original_data['data'] = -original_data.iloc[:,2]
    original_data.set_index('orf', inplace=True)
    original_data = original_data[['data']].copy()   
    original_data = original_data.groupby(original_data.index).mean()
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[8]:


original_data = pd.concat(original_data_list, axis=0)


# In[10]:


original_data.shape


# In[11]:


original_data = original_data.groupby(original_data.index).mean()


# In[12]:


original_data.shape


# # Load & process tested strains

# In[13]:


tested = pd.read_excel('raw_data/CompleteDeletionLibrary.xlsx', sheet_name='Sheet1')


# In[14]:


tested.head()


# In[15]:


tested['orf'] = tested['Unique Name'].astype(str)


# In[16]:


tested['orf'] = clean_orf(tested['orf'])


# In[19]:


typo_fixes = {'YAR002AW':'YAR002W-A','YOLO57W':'YOL057W','YKLO72W':'YKL072W','YJL206-A':'YJL206C-A',
              'YLR287-A':'YLR287C-A','YFL033AC':'YFL033C-A','YOLO62C':'YOL062C'}


# In[21]:


tested['orf'] = tested['orf'].apply(lambda x: typo_fixes[x] if x in typo_fixes.keys() else x)


# In[22]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[23]:


t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[24]:


tested_orfs = np.unique(tested['orf'].values)


# In[25]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[26]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[27]:


data = original_data.copy()


# In[28]:


dataset_ids = [569]
datasets = datasets.reindex(index=dataset_ids)


# In[29]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[30]:


data.head()


# ## Subset to the genes currently in SGD

# In[31]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[32]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[33]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[34]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[35]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

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




