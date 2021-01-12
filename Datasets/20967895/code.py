#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 20967895
paper_name = 'yadav_yadav_2011' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/Workbook1yea_1825_supportinginfor.xlsx', sheet_name='Sheet1')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[9]:


original_data_list = []
for c in original_data.columns:
    genes = original_data.loc[0,c].split(',')
    dt = pd.DataFrame(data={'gene': genes})
    dt['gene'] = dt['gene'].astype(str)
    dt['gene'] = clean_genename(dt['gene'])
    
    dt.loc[dt['gene']=='FMP53','gene'] = 'COQ9'
    dt.loc[dt['gene']=='PR116W','gene'] = 'YPR116W'
    dt.loc[dt['gene']=='GON5','gene'] = 'YPL183W-A'
    
    dt['orf'] = translate_sc(dt['gene'], to='orf')
    t = looks_like_orf(dt['orf'])
    print(dt.loc[~t,])
    
    dt = dt.loc[t,:]
    dt['data'] = -1
    dt.set_index('orf', inplace=True)
    dt = dt.groupby(dt.index).mean()
    
    original_data_list.append(dt)


# In[10]:


original_data = pd.concat(original_data_list, axis=1)


# In[11]:


original_data.head()


# In[12]:


original_data.index.name = 'orf'


# # Prepare the final dataset

# In[13]:


data = original_data.copy()


# In[14]:


dataset_ids = [132,446,447]
datasets = datasets.reindex(index=dataset_ids)


# In[15]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[16]:


data.head()


# ## Subset to the genes currently in SGD

# In[17]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[18]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[19]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[20]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[21]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[22]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[23]:


from IO.save_data_to_db3 import *


# In[24]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




