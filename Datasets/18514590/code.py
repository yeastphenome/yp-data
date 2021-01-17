#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 18514590
paper_name = 'serero_boiteux_2008' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_csv('raw_data/hits_genenames.txt', header=None, sep='\t')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


original_data.head()


# In[8]:


original_data['gene'] = original_data[0].astype(str)


# In[9]:


# Eliminate all white spaces & capitalize
original_data['gene'] = clean_orf(original_data['gene'])


# In[10]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['gene'], to='orf')


# In[11]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[13]:


original_data.tail()


# In[15]:


original_data['data'] = original_data[1].apply(lambda x: -(len(x)+1) if isinstance(x, str) else -1)


# In[16]:


original_data.set_index('orf', inplace=True)


# In[17]:


original_data = original_data[['data']].copy()


# In[18]:


original_data = original_data.groupby(original_data.index).mean()


# In[19]:


original_data.shape


# # Load & process tested strains

# In[ ]:


# Attempt 1: go through the EXCLE files. PROBLEM: too many genes obtained (5617)


# In[ ]:


# Attempt 2: go through the DOC files. Covert all DOC files into TXT files by running: sudo textutil -convert txt */*.DOC. 
# Read the TXT files that end with 1


# In[76]:


txt_files = [f for f in os.listdir('raw_data/') if os.path.isfile(os.path.join('raw_data/', f)) 
             and (f.endswith('~1.txt') or f.endswith('a.txt'))]


# In[77]:


len(txt_files)


# In[78]:


tested_orfs = []
for f in txt_files:
    t = pd.read_csv('raw_data/' + f, header=None, sep='\t')
    tested_orfs.append(t)


# In[79]:


tested = pd.concat(tested_orfs, axis=0, ignore_index=True)


# In[80]:


tested['orf'] = tested[0].astype(str)


# In[81]:


tested['orf'] = clean_orf(tested['orf'])


# In[82]:


tested['orf'] = translate_sc(tested['orf'].values, to='orf')


# In[83]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[84]:


tested = tested.loc[t,:]


# In[85]:


tested_orfs = tested['orf'].unique()


# In[86]:


tested_orfs.shape


# In[88]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[89]:


len(missing)


# In[90]:


tested_orfs = list(tested_orfs) + missing


# In[91]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[92]:


data = original_data.copy()


# In[93]:


dataset_ids = [99]
datasets = datasets.reindex(index=dataset_ids)


# In[94]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[95]:


data.head()


# ## Subset to the genes currently in SGD

# In[96]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[97]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[98]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[99]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[100]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[101]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[102]:


from IO.save_data_to_db3 import *


# In[103]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




