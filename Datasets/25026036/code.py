#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 25026036
paper_name = 'elbaz_alon_schuldiner_2014' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[8]:


original_data = pd.read_excel('raw_data/mmc2.xlsx', sheet_name='Table S1', skiprows=1)


# In[9]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[10]:


original_data.head()


# In[11]:


original_data['orf'] = original_data['ORF'].astype(str)


# In[12]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[13]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[14]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[16]:


phenotypes = original_data['Phenotype'].unique()
phenotypes


# In[17]:


phenotypes_digit = dict()
phenotypes_digit['Decrease in GFP intensity compared to wt'] = -1
phenotypes_digit['GFP foci less defined compard to wt'] = -1
phenotypes_digit['Increase in GFP intensity compared to wt'] = 1
phenotypes_digit['Increase in amount of foci compared to wt'] = 1
phenotypes_digit['No MCS, signal spread al over mitochondria'] = -2


# In[18]:


original_data['data'] = original_data['Phenotype'].apply(lambda x: phenotypes_digit[x])


# In[19]:


original_data.set_index('orf', inplace=True)


# In[20]:


original_data = original_data[['data']].copy()


# In[21]:


original_data = original_data.groupby(original_data.index).mean()


# In[22]:


original_data.shape


# # Load & process tested strains

# In[23]:


tested = pd.read_excel('raw_data/KO_DAmP_ORFs.xlsx', sheet_name='Sheet1', skiprows=1)


# In[24]:


tested.head()


# In[26]:


tested['orf'] = tested['ORF '].astype(str)


# In[27]:


tested['orf'] = clean_orf(tested['orf'])


# In[30]:


typo_fixes = {'YOLO57W':'YOL057W','YOLO62C':'YOL062C','YBRF182C-A':'YBR182C-A','YLR287-A':'YLR287C-A','YJL206-A':'YJL206C-A'}


# In[31]:


tested['orf'] = tested['orf'].apply(lambda x: typo_fixes[x] if x in typo_fixes.keys() else x)


# In[32]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[33]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[34]:


tested = tested.loc[t,:]


# In[35]:


tested_orfs = np.unique(tested['orf'].values)


# In[36]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[37]:


# Removing the missing strains (they were tested as DAMP strains, not deletions?)
original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# In[38]:


original_data.head()


# # Prepare the final dataset

# In[39]:


data = original_data.copy()


# In[40]:


dataset_ids = [16000]
datasets = datasets.reindex(index=dataset_ids)


# In[41]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[42]:


data.head()


# ## Subset to the genes currently in SGD

# In[43]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[44]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[45]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[46]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[47]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[48]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[49]:


from IO.save_data_to_db3 import *


# In[50]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:



