#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 18688276
paper_name = 'ericson_nislow_2008' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[109]:


original_data = pd.read_excel('raw_data/pgen.1000151.s003.xlsx', sheet_name='Sheet1', skiprows=2)


# In[110]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[111]:


original_data.head()


# In[112]:


original_data['orf'] = original_data['ORF'].astype(str)


# In[113]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[114]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[115]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[116]:


original_data.set_index('orf', inplace=True)


# In[117]:


data_cols = original_data.columns[np.arange(2,88)]


# In[118]:


original_data1 = original_data.loc[original_data['Zygosity']=='hom',data_cols].copy()
original_data2 = original_data.loc[original_data['Zygosity']=='het',data_cols].copy()


# In[119]:


original_data1 = original_data1.apply(pd.to_numeric, axis=1, errors='coerce')
original_data2 = original_data2.apply(pd.to_numeric, axis=1, errors='coerce')


# In[120]:


original_data1 = original_data1.groupby(original_data1.index).mean()
original_data2 = original_data2.groupby(original_data2.index).mean()


# In[121]:


original_data1.shape


# In[122]:


original_data2.shape


# # Load dataset info

# In[123]:


original_data1.columns = [c.split('.')[0] for c in original_data1.columns.values]
original_data2.columns = [c.split('.')[0] for c in original_data2.columns.values]


# In[124]:


df = pd.read_excel('extras/drug_dataset.xlsx', sheet_name='Sheet1')
df.set_index('Unnamed: 0', inplace=True)


# In[125]:


df.drop_duplicates(subset=['Hom dataset','Het dataset'], inplace=True)
df.shape


# In[126]:


hom_dataset_id = [df.loc[c,'Hom dataset'] for c in original_data1.columns.values]
het_dataset_id = [df.loc[c,'Het dataset'] for c in original_data1.columns.values]


# In[127]:


original_data1.columns = hom_dataset_id
original_data2.columns = het_dataset_id


# In[128]:


original_data1 = original_data1.T
original_data1 = original_data1.groupby(original_data1.index).mean()
original_data1 = original_data1.T
original_data1.shape


# In[129]:


original_data2 = original_data2.T
original_data2 = original_data2.groupby(original_data2.index).mean()
original_data2 = original_data2.T
original_data2.shape


# In[130]:


original_data = original_data1.join(original_data2, how='outer')


# In[131]:


original_data.head()


# In[132]:


# Taking the opposite because the original data is log2(ctrl/treatment)
original_data = -original_data


# # Prepare the final dataset

# In[133]:


data = original_data.copy()


# In[134]:


dataset_ids = original_data.columns.values
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




