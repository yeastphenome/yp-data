#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 15525520
paper_name = 'pan_boeke_2004' 


# In[91]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[92]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[97]:


original_data = pd.read_excel('raw_data/TableS2.xlsx', sheet_name='Table 1', skiprows=1)


# In[98]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[99]:


original_data.head()


# In[100]:


original_data['gene'] = original_data['Gene Name'].astype(str)


# In[101]:


# Eliminate all white spaces & capitalize
original_data['gene'] = clean_genename(original_data['gene'])


# In[102]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['gene'], to='orf')


# In[103]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[104]:


original_data = original_data.loc[t,:]


# In[105]:


original_data.set_index('orf', inplace=True)


# In[106]:


original_data.drop(columns=['Gene Name','gene'], inplace=True)


# In[107]:


# Reversing the sign to indicate that lower numbers = slower growth
original_data = -original_data.apply(pd.to_numeric, axis=1, errors='coerce')


# In[108]:


original_data = original_data.groupby(original_data.index).mean()


# In[109]:


original_data.shape


# # Load & process tested strains

# In[110]:


# Microarray platform
tested = pd.read_csv('raw_data/GPL1444.txt', sep='\t', skiprows=273)


# In[111]:


tested.head()


# In[112]:


# Untreated control
tested2 = pd.read_csv('raw_data/GSM30549.txt', sep='\t', skiprows=52)


# In[113]:


tested2.head()


# In[114]:


# Get all ORFS that have a flag == 0 (passed the basic filter)
tested = tested.merge(tested2[['ID_REF','Flags']], left_on='ID', right_on='ID_REF', how='left')


# In[115]:


tested = tested.loc[tested['Flags'] > -50,:]


# In[116]:


tested['orf'] = tested['ORF'].astype(str)


# In[117]:


tested['orf'] = clean_orf(tested['orf'])


# In[118]:


tested['orf'] = translate_sc(tested['orf'].values, to='orf')


# In[119]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[120]:


tested = tested.loc[t,:]


# In[121]:


tested_orfs = tested['orf'].unique()


# In[122]:


# Excluded strains
tested3 = pd.read_excel('raw_Data/TableS3.xlsx', sheet_name='Table 1', skiprows=1)
tested3.head()


# In[123]:


tested3['orf'] = tested3['ORF name'].astype(str)
tested3['orf'] = clean_orf(tested3['orf'])
tested3['orf'] = translate_sc(tested3['orf'].values, to='orf')


# In[124]:


# Make sure everything translated ok
t = looks_like_orf(tested3['orf'])
print(tested3.loc[~t,])


# In[125]:


tested3 = tested3.loc[t,:]


# In[126]:


excluded = tested3['orf'].values


# In[127]:


tested_orfs = [orf for orf in tested_orfs if orf not in excluded]


# In[128]:


len(excluded)


# In[129]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[130]:


tested_orfs = tested_orfs + missing


# In[131]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[132]:


data = original_data.copy()


# In[133]:


dataset_ids = np.arange(5237, 5246)
datasets = datasets.reindex(index=dataset_ids)


# In[135]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[136]:


data.head()


# ## Subset to the genes currently in SGD

# In[137]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[138]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[139]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[140]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[141]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[142]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[143]:


from IO.save_data_to_db3 import *


# In[144]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




