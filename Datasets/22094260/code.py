#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 22094260
paper_name = 'skrtic_schimmer_2011' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[27]:


original_data = pd.read_csv('raw_data/het_damp.rawsummary', sep='\t')


# In[28]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[29]:


original_data.head()


# In[30]:


# First, eliminate the data for the DAMP strains
original_data = original_data.loc[~original_data['Hybridization REF'].str.contains('DAMP'),]
original_data.shape


# In[31]:


# Now, extract the ORF
original_data['orf'] = original_data['Hybridization REF'].apply(lambda x: x[0:x.find(':')])


# In[32]:


original_data['orf'] = original_data['orf'].astype(str)


# In[33]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[34]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[35]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[36]:


original_data = original_data.loc[t,]


# In[37]:


original_data.set_index('orf', inplace=True)


# In[38]:


for c in original_data.columns.values[1:]:
    original_data[c] = pd.to_numeric(original_data[c])


# In[39]:


# Take the average of the 2 YPGE_DMSO controls
original_data['YPGE_DMSO_avg'] = original_data[['10_11_24_YPGE_DMSO','10_11_24_YPGE_DMSO_2']].mean(axis=1)


# In[40]:


# Divide each treatment by its control
original_data['10_11_24_YPGE_chloramph_0.79_norm'] = original_data['10_11_24_YPGE_chloramph_0.79'] / original_data['YPGE_DMSO_avg']
original_data['10_11_24_YPGE_chloramph_0.99_norm'] = original_data['10_11_24_YPGE_chloramph_0.99'] / original_data['YPGE_DMSO_avg']
original_data['10_11_24_YPGE_doxorub.12.5_norm'] = original_data['10_11_24_YPGE_doxorub.12.5'] / original_data['YPGE_DMSO_avg']
original_data['10_11_24_YPGE_linezol_47.1_norm'] = original_data['10_11_24_YPGE_linezol_47.1'] / original_data['YPGE_DMSO_avg']

original_data['10_12_10_tigecyc51.5uM_norm'] = original_data['10_12_10_tigecyc51.5uM'] / original_data['10_12_10_tigecycDMSOctrl']
original_data['10_12_10_tigecyc64.4uM_norm'] = original_data['10_12_10_tigecyc64.4uM'] / original_data['10_12_10_tigecycDMSOctrl']
original_data['10_12_10_tigecyc80.5uM_norm'] = original_data['10_12_10_tigecyc80.5uM'] / original_data['10_12_10_tigecycDMSOctrl']


# In[41]:


cols_to_keep = ['10_11_24_YPGE_chloramph_0.79_norm','10_11_24_YPGE_chloramph_0.99_norm',
                '10_11_24_YPGE_doxorub.12.5_norm',
                '10_11_24_YPGE_linezol_47.1_norm',
                '10_12_10_tigecyc51.5uM_norm','10_12_10_tigecyc64.4uM_norm','10_12_10_tigecyc80.5uM_norm']


# In[42]:


original_data = original_data[cols_to_keep]


# In[43]:


original_data = original_data.groupby(original_data.index).mean()


# In[44]:


original_data.shape


# In[45]:


original_data.head()


# # Prepare the final dataset

# In[46]:


data = original_data.copy()


# In[47]:


dataset_ids = [16572,16591,16570,16573,16571,16592,16593]
datasets = datasets.reindex(index=dataset_ids)


# In[48]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[49]:


data.head()


# ## Subset to the genes currently in SGD

# In[50]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[51]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[52]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[53]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[54]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[55]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[56]:


from IO.save_data_to_db3 import *


# In[57]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




