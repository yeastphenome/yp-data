#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[46]:


paper_pmid = 32064787
paper_name = 'mattiazzi_usaj_andrews_2020' 


# In[47]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[48]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# In[96]:


original_data = pd.read_excel('raw_data/msb199243-sup-0003-tableev2.xlsx', sheet_name='penetrance_phenotype_data')


# In[97]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[98]:


original_data['ORF'] = original_data['ORF'].astype(str)


# In[99]:


# Eliminate all white spaces & capitalize
original_data['ORF'] = clean_orf(original_data['ORF'])


# In[100]:


original_data['ORF'] = translate_sc(original_data['ORF'], to='orf')


# In[101]:


# Make sure everything translated ok
t = looks_like_orf(original_data['ORF'])


# In[102]:


print(original_data.loc[~t,])


# In[103]:


# Eliminate strains that are not deletions
dels = original_data['StrainID'].str.startswith('DMA')
original_data = original_data.loc[dels.values,:]
print(original_data.shape)


# In[104]:


original_data.set_index('ORF', inplace=True)
original_data.index.name='orf'


# In[105]:


dataset_map = {'actin_aggregate': 16403,
               'actin_bright_patches': 16415,
               'actin_decreased_patch_number': 16416,
               'actin_depolarized_patches': 16417,
               'coat_aggregate': 16418,
               'coat_decreased_patch_number': 16419,
               'coat_depolarized_patches': 16420,
               'coat_increased_patch_number': 16421,
               'LE_fragmented': 16422,
               'LE_membrane': 16423,
               'LE_condensed': 16424,
               'vacuole_class_E': 16425,
               'vacuole_enlarged': 16426,
               'vacuole_VATPase_defect': 16427,
               'vacuole_fragmented': 16428,
               'vacuole_class_G': 16429,
               'vacuole_multilobed': 16430}


# In[106]:


# Get the relevant columns
original_data = original_data.loc[:, dataset_map.keys()]


# In[107]:


print(original_data.shape)


# In[108]:


# If the same strain is present more than once, average its values
original_data = original_data.groupby(original_data.index).mean()


# In[109]:


print('Final data dimensions: %d x %d' % (original_data.shape))


# # Prepare the final dataset

# In[111]:


data = original_data.copy()


# In[112]:


dataset_ids = [dataset_map[c] for c in data.columns.values]


# In[113]:


datasets = datasets.reindex(index=dataset_ids)


# In[114]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# ## Subset to the genes currently in SGD

# In[116]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[118]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# # Normalize

# In[119]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[120]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[121]:


data_norm[data.isnull()] = np.nan


# In[122]:


data_all = data.join(data_norm)


# # Print out

# In[123]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[124]:


from IO.save_data_to_db3 import *


# In[125]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




