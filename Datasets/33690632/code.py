#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 33690632
paper_name = 'nicastro_devirgilio_2021' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


file = 'raw_data/Rapa-IAA_full.xlsx'
sheets = ['2017-04-14-IAA','2017-04-03-Rapamycin']


# In[6]:


original_data_list = []
for s in sheets:
    original_data = pd.read_excel(file, sheet_name=s)
    print('Original data dimensions: %d x %d' % (original_data.shape))
    original_data['orf'] = original_data['ID Column'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data.set_index('orf', inplace=True)
    original_data['data'] = original_data['Growth Ratio (Comparer / Exp)'].astype(float)
    original_data['data'] = 1/original_data['data']
    original_data = original_data[['data']].copy()
    original_data = original_data.groupby(original_data.index).mean()
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[7]:


original_data = pd.concat(original_data_list, axis=1)


# In[8]:


original_data.index.name = 'orf'


# In[9]:


original_data.head()


# In[10]:


original_data.shape


# In[12]:


# Remove DAMP strains from the list
damp = pd.read_excel('raw_data/YSC5090 - Yeast DAmP HAPLOID_v2.xls', sheet_name='data')


# In[14]:


damp['orf'] = damp['ORF'].astype(str)
damp['orf'] = clean_orf(damp['orf'])
damp['orf'] = translate_sc(damp['orf'], to='orf')


# In[16]:


t = looks_like_orf(damp['orf'])
print(damp.loc[~t,])


# In[17]:


damp_orfs = damp['orf'].unique()


# In[21]:


damp_orfs = [orf for orf in damp_orfs if orf in original_data.index.values]
len(damp_orfs)


# In[22]:


original_data = original_data.drop(index=damp_orfs)


# In[23]:


original_data.shape


# # Prepare the final dataset

# In[24]:


data = original_data.copy()


# In[25]:


dataset_ids = [22056,22055]
datasets = datasets.reindex(index=dataset_ids)


# In[26]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[27]:


data.head()


# ## Subset to the genes currently in SGD

# In[28]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[29]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[30]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[31]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[32]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[33]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[34]:


from IO.save_data_to_db3 import *


# In[35]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:



