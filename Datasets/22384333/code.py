#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 22384333
paper_name = 'hoon_nislow_2011' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[6]:


files = ['TableS1.xlsx','TableS2.xlsx','TableS3.xlsx']


# In[12]:


original_data_list = []
for f in files:
    xl = pd.ExcelFile('raw_data/' + f)
    sheets = xl.sheet_names
    original_data = pd.read_excel('raw_data/' + f, sheet_name=sheets[0], skiprows=3)
    print('Original data dimensions: %d x %d' % (original_data.shape))
#     print(original_data.head())
    
    if 'ORF' in original_data.columns:
        c = 'ORF'
        data_col = 'Log2 ratio'
    else:
        c = 'strain'
        data_col = 'Average'
    
    original_data['orf'] = original_data[c].apply(lambda x: x.split(':')[0])
    
    # Eliminate all white spaces & capitalize
    original_data['orf'] = clean_orf(original_data['orf'])
    
    # Translate to ORFs 
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    
    # Make sure everything translated ok
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    
    original_data['data'] = original_data[data_col]
    if c == 'strain':
        original_data['data'] = -original_data['data']
        
    original_data.set_index('orf', inplace=True)
    original_data = original_data[['data']].copy()
    original_data = original_data.groupby(original_data.index).mean()
    
    print(original_data.shape)
    original_data_list.append(original_data)


# In[13]:


original_data = pd.concat(original_data_list, axis=1)


# In[15]:


original_data.index.name='orf'


# In[20]:


dt = original_data.iloc[:,[0,1]].values
dt[np.isnan(dt)] = 0
original_data.iloc[:,[0,1]] = dt


# In[21]:


original_data.head()


# # Prepare the final dataset

# In[22]:


data = original_data.copy()


# In[23]:


dataset_ids = [519, 755, 16253]
datasets = datasets.reindex(index=dataset_ids)


# In[24]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[25]:


data.head()


# ## Subset to the genes currently in SGD

# In[26]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[27]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[28]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[29]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[30]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[31]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[32]:


from IO.save_data_to_db3 import *


# In[33]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




