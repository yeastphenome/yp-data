#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 11907266
paper_name = 'dimmer_westermann_2002' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[18]:


files = ['no_growth.txt','poor_growth.txt']
scores = [-1, -0.5]


# In[19]:


original_data_list = []
for ixf, f in enumerate(files):
    original_data = pd.read_csv('raw_data/' + f, sep='\n', header=None)
    print('Original data dimensions: %d x %d' % (original_data.shape))
    print(original_data.head())
    
    original_data['orf'] = original_data[0].apply(lambda x: re.split(' |\t',x)[0])
    original_data['orf'] = original_data['orf'].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    
    original_data = original_data.loc[t,:]
    original_data['data'] = scores[ixf]
    original_data.set_index('orf', inplace=True)
    original_data = original_data[['data']].copy()
    original_data = original_data.groupby(original_data.index).mean()
    
    print(original_data.shape)
    
    original_data_list.append(original_data)


# In[22]:


original_data = pd.concat(original_data_list, axis=0)
original_data = original_data.groupby(original_data.index).mean()


# In[23]:


original_data.head()


# In[24]:


original_data.shape


# # Load & process tested strains

# In[55]:


tested = pd.read_csv('raw_data/HOMOZYGOUS DIPLOID 1+2 ResGen.txt', sep=' ')


# In[56]:


tested.head()


# In[57]:


tested['orf'] = tested['Unnamed: 3'].astype(str)


# In[58]:


tested['orf'] = clean_orf(tested['orf'])


# In[59]:


tested.loc[tested['orf']=='YELOO1C','orf'] = 'YEL001C'


# In[60]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[61]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[62]:


tested = tested.loc[t,:]


# In[63]:


tested_orfs = tested['orf'].unique()


# In[64]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[65]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[66]:


data = original_data.copy()


# In[67]:


dataset_ids = [470]
datasets = datasets.reindex(index=dataset_ids)


# In[68]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[69]:


data.head()


# ## Subset to the genes currently in SGD

# In[70]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[71]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[72]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[73]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[74]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[75]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[76]:


from IO.save_data_to_db3 import *


# In[77]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:



