#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 17846143
paper_name = 'morton_coote_2007' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/Table1.xlsx', sheet_name='Sheet1')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[9]:


d_genes = []
mag_genes = []

cols = original_data.columns.values[2:]
for r in original_data.iterrows():
    for c in cols:
        s = str(r[1][c])
        s = s.replace('\xa0','')
        genes = s.split(',')
        if not isinstance(genes, list):
            genes = [genes]
        if c == cols[0]:
            d_genes = d_genes + genes
            mag_genes = mag_genes + genes
        elif c == cols[1]:
            d_genes = d_genes + genes
        elif c == cols[2]:
            mag_genes = mag_genes + genes


# In[10]:


d_genes = [s.strip() for s in d_genes if not s == 'nan']
mag_genes = [s.strip() for s in mag_genes if not s == 'nan']


# In[11]:


d_genes = clean_genename(d_genes)
mag_genes = clean_genename(mag_genes)


# In[12]:


d_orfs = translate_sc(d_genes, to='orf')
mag_orfs = translate_sc(mag_genes, to='orf')


# In[13]:


d_orfs = np.array(d_orfs)
mag_orfs = np.array(mag_orfs)


# In[14]:


d_orfs[d_orfs=='TMA29'] = 'YMR226C'


# In[15]:


t = looks_like_orf(d_orfs)
print(d_orfs[~np.array(t)])


# In[16]:


t = looks_like_orf(mag_orfs)
print(mag_orfs[~np.array(t)])


# In[17]:


all_orfs = np.unique(np.concatenate((d_orfs, mag_orfs)))


# In[18]:


data = pd.DataFrame(index=all_orfs, columns=['D','M'], data=np.zeros((len(all_orfs),2)))


# In[19]:


data.loc[d_orfs,'D'] = -1
data.loc[mag_orfs,'M'] = -1


# In[21]:


data.index.name = 'orf'


# # Prepare the final dataset

# In[22]:


dataset_ids = [16536,16535]
datasets = datasets.reindex(index=dataset_ids)


# In[23]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[24]:


data.head()


# In[25]:


data = data.groupby(data.index).mean()


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

# In[29]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[30]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[31]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[32]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[33]:


from IO.save_data_to_db3 import *


# In[34]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




