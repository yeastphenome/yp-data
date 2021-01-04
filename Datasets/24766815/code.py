#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 24766815
paper_name = 'kemmeren_holstege_2014' 


# In[125]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[126]:


datasets.set_index('dataset_id', inplace=True)


# In[128]:


datasets['orf'] = datasets['name'].apply(lambda x: x[x.find("(")+1:x.find(")")])


# In[131]:


datasets.shape


# In[129]:


datasets.head()


# # Load & process the data

# In[97]:


original_data = pd.read_csv('raw_data/deleteome_all_mutants_ex_wt_var_controls.txt', sep='\t', low_memory=False)


# In[98]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[99]:


original_data = original_data.T


# In[100]:


original_data.columns = original_data.loc['systematicName',:].astype(str)


# In[101]:


original_data.drop(index=['reporterId','systematicName','geneSymbol'], inplace=True)


# In[102]:


data_rows = original_data.index[original_data.iloc[:,0]=='M']


# In[103]:


original_data = original_data.loc[data_rows,:]


# In[104]:


original_data.drop(columns=['nan'], inplace=True)


# In[105]:


original_data.head()


# In[106]:


# Extract deletion strain names
original_data['genes'] = [x.split('-del')[0] for x in original_data.index.values]


# In[107]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[108]:


manual_fixes = {'LUG1': 'YLR352W','CYCC': 'YNL025C'}
original_data['genes'] = original_data['genes'].apply(lambda x: manual_fixes[x] if x in manual_fixes.keys() else x)


# In[109]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['genes'], to='orf')


# In[110]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[111]:


original_data = original_data.loc[t,:]


# In[112]:


original_data.set_index('orf', inplace=True)
original_data.drop(columns=['genes'], inplace=True)


# In[113]:


original_data = original_data.astype(float)


# In[114]:


original_data = original_data.groupby(original_data.index).mean()


# In[115]:


original_data = original_data.T
original_data = original_data.groupby(original_data.index).mean()
original_data = original_data.T


# In[116]:


original_data.shape


# In[117]:


original_data.head()


# # Prepare the final dataset

# In[141]:


data = original_data.copy()


# In[142]:


datasets = datasets.reset_index().set_index('orf')
dataset_ids = datasets.reindex(index=original_data.columns.values)['dataset_id'].values


# In[143]:


datasets = datasets.reset_index().set_index('dataset_id')


# In[144]:


datasets = datasets.reindex(index=dataset_ids)


# In[145]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[146]:


data.head()


# ## Subset to the genes currently in SGD

# In[147]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[148]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[152]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[153]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[154]:


data_vals = data.values
data_norm_vals = data_norm.values
data_norm_vals[np.isnan(data_vals)] = np.nan


# In[155]:


data_norm = pd.DataFrame(index=data_norm.index, columns=data_norm.columns, data=data_norm_vals)


# In[156]:


data_all = data.join(data_norm)
data_all.head()


# # Print out

# In[157]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[162]:


# from IO.save_data_to_db3 import *


# In[163]:


# save_data_to_db(data_all, paper_pmid)


# In[170]:


data_long = pd.melt(data.droplevel('orf', axis=0).reset_index(), id_vars=['gene_id'], col_level='dataset_id')


# In[172]:


data_norm_long = pd.melt(data_norm.droplevel('orf', axis=0).reset_index(), id_vars=['gene_id'], col_level='dataset_id')


# In[174]:


data_long = data_long.merge(data_norm_long, on=['gene_id','dataset_id'])


# In[176]:


data_long.columns = ['gene_id','dataset_id','value','valuez']


# In[177]:


data_long.to_csv('kemmeren_holstege_2014_long.txt', sep='\t', index=False)


# In[ ]:




