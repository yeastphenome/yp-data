#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 12496357
paper_name = 'begley_samson_2002' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[72]:


original_data = pd.read_excel('raw_data/ORIG130404_Begley2001raw.xlsx', sheet_name='Sheet1')


# In[73]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[74]:


original_data.head()


# In[75]:


original_data['orf'] = original_data['ORF'].astype(str)


# In[76]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[77]:


typo_fixes = {'YKLO72W':'YKL072W','YOLO57W':'YOL057W','YOLO62C':'YOL062C'}
original_data['orf'] = original_data['orf'].apply(lambda x: typo_fixes[x] if x in typo_fixes.keys() else x)


# In[78]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[79]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[80]:


original_data = original_data.loc[t,:]


# In[81]:


original_data['Treatment'].unique()


# In[82]:


data_cols = [c for c in original_data.columns if 'control/dose' in c]
data_cols


# In[83]:


original_data.set_index('orf', inplace=True)


# In[84]:


original_data = original_data[['Treatment'] + data_cols].copy()


# In[85]:


original_data[data_cols] = original_data[data_cols].apply(pd.to_numeric, axis=1, errors='coerce')


# In[86]:


original_data[data_cols] = 1 / original_data[data_cols]


# In[87]:


original_data.head()


# In[88]:


for dose in np.arange(1,5):
    cols = [c for c in original_data.columns if 'dose ' + str(dose) in c]
    original_data['dose' + str(dose)] = original_data[cols].mean(axis=1)


# In[89]:


treatments = original_data['Treatment'].unique()


# In[90]:


original_data_list = []
for t in treatments:
    original_data1 = original_data.loc[original_data['Treatment']==t,['dose1','dose2','dose3','dose4']].copy()
    original_data1 = original_data1.groupby(original_data1.index).mean()
    cols = [t+'_'+c for c in original_data1.columns]
    original_data1.columns=cols
    
    original_data_list.append(original_data1)


# In[91]:


original_data = pd.concat(original_data_list, axis=1)


# In[92]:


original_data.shape


# In[93]:


original_data.head()


# In[94]:


dt = pd.read_csv('extras/phenotype_datasetids.txt', sep='\t', header=None)


# In[96]:


dataset_ids = dt[1].values


# In[97]:


original_data.shape


# # Prepare the final dataset

# In[98]:


data = original_data.copy()


# In[99]:


datasets = datasets.reindex(index=dataset_ids)


# In[100]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[101]:


data.head()


# ## Subset to the genes currently in SGD

# In[102]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[103]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[104]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[105]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[106]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[107]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[108]:


from IO.save_data_to_db3 import *


# In[109]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




