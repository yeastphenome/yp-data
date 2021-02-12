#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 12871900
paper_name = 'griffith_devine_2003' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[21]:


original_data = pd.read_csv('raw_data/hits.txt', header=None, sep='\t')


# In[22]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[23]:


original_data.head()


# In[24]:


original_data['gene'] = original_data[0].astype(str)


# In[25]:


original_data['gene'] = original_data['gene'].apply(lambda x: x.split(' ')[0])


# In[26]:


# Eliminate all white spaces & capitalize
original_data['gene'] = clean_genename(original_data['gene'])


# In[27]:


fixes = {'TCI1': 'YDR161W', 'SDF1': 'YPR040W'}
original_data['gene'] = original_data['gene'].apply(lambda x: fixes[x] if x in fixes.keys() else x)


# In[28]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['gene'], to='orf')


# In[29]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[30]:


original_data = original_data.loc[t,:]


# In[31]:


original_data['data'] = pd.to_numeric(original_data[1], errors='coerce')


# In[32]:


original_data.set_index('orf', inplace=True)


# In[33]:


original_data = original_data[['data']].copy()


# In[34]:


original_data = original_data.groupby(original_data.index).mean()


# In[35]:


original_data.shape


# # Load & process tested strains

# In[37]:


files = ['Res Gen diploid knock01.xlsx','Res Gen diploid knockouts02.xlsx']
sheets = ['Res Gen diploid knock01.txt','Res Gen diploid knockouts02.txt']


# In[38]:


tested_list = []
for ixf, f in enumerate(files):
    tested = pd.read_excel('raw_data/' + f, sheet_name=sheets[ixf], skiprows=1)
#     print(tested.head())
    tested['orf'] = tested['ORF name'].astype(str)
    tested['orf'] = clean_orf(tested['orf'])
    
    typo_fixes = {'TAL004W':'YAL004W','YELOO1C':'YEL001C','KL187C':'YKL187C'}
    tested['orf'] = tested['orf'].apply(lambda x: typo_fixes[x] if x in typo_fixes.keys() else x)
    tested['orf'] = translate_sc(tested['orf'], to='orf')
    # Make sure everything translated ok
    t = looks_like_orf(tested['orf'])
    print(tested.loc[~t,])
    tested = tested.loc[t,:]
    
    if 'plate' in tested.columns:
        c = 'plate'
    else:
        c = 'row'
    tested = tested.loc[(tested[c] >= 301) & (tested[c] <= 349)]
    tested_list.append(tested[['orf']])


# In[39]:


tested = pd.concat(tested_list, axis=0)


# In[40]:


tested.head()


# In[41]:


tested_orfs = tested['orf'].unique()


# In[42]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[43]:


tested_orfs = list(tested_orfs) + missing


# In[44]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[45]:


data = original_data.copy()


# In[46]:


dataset_ids = [480]
datasets = datasets.reindex(index=dataset_ids)


# In[47]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[48]:


data.head()


# ## Subset to the genes currently in SGD

# In[49]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[50]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[51]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[52]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[53]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[54]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[55]:


from IO.save_data_to_db3 import *


# In[56]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




