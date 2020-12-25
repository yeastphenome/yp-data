#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 29507053
paper_name = 'salignon_yvert_2018' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data (part 1)

# In[5]:


# Data obtained from R by running:
# > load(file="sup_file_Met.rda")
# > write.table(tfit.summary, file="Met_tfit_summary.txt", sep='\t', quote=FALSE, row.names=FALSE)


# In[6]:


# If ORF ends with a "1", strip the "1"
def remove_trailing_one(s):
    return s[:-1] if s[-1]=='1' else s


# In[7]:


original_data_list = []
files = ['Salt_tfit_summary.txt','Met_tfit_summary.txt']


# In[8]:


original_data_list2 = []
for f in files:
    original_data = pd.read_csv('raw_data/' + f, sep='\t')
    print('Original data dimensions: %d x %d' % (original_data.shape))
    original_data['orf'] = original_data['orf'].astype(str)
    
    # Eliminate all white spaces & capitalize
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = original_data['orf'].apply(remove_trailing_one)
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    
    cols = ['orf','w.N','w.S','w.6h','w.12h','w.18h','w.24h','w.42h']
    original_data = original_data.loc[:,cols]
    original_data.set_index('orf', inplace=True)
    original_data = original_data.groupby(original_data.index).mean()
    print(original_data.shape)
    
    original_data_list2.append(original_data)


# In[9]:


original_data1, original_data2 = original_data_list2


# In[10]:


original_data = original_data1.join(original_data2, how='outer', lsuffix="_1", rsuffix="_2")


# In[11]:


original_data.head()


# In[12]:


dataset_ids = [16167, 16168, 16169, 16175, 16176, 16177, 16178] + [16170, 16171, 16172, 16179, 16180, 16181, 16182]


# # Prepare the final dataset

# In[13]:


data = original_data.copy()


# In[14]:


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


# In[19]:


data.head()


# # Normalize

# In[20]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[21]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[22]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[23]:


data_all.head()


# # Print out

# In[24]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[25]:


from IO.save_data_to_db3 import *


# In[26]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




