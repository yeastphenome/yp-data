#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 16467474
paper_name = 'chasse_dohlman_2006' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[76]:


original_data = pd.read_excel('raw_data/Table_S1_copy.xlsx', sheet_name='Table 1')


# In[77]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[78]:


original_data.head()


# In[79]:


original_data1 = original_data.iloc[:,[0,1,2]].copy()
original_data2 = original_data.iloc[:,[3,4,5]].copy()


# In[80]:


original_data2.columns = original_data1.columns


# In[81]:


original_data = pd.concat([original_data1, original_data2], axis=0)


# In[82]:


original_data['orf'] = original_data['Gene'].astype(str)


# In[83]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[84]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'].values, to='orf')


# In[85]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[86]:


original_data = original_data.loc[t,:]


# In[87]:


original_data = original_data.reset_index()


# In[88]:


original_data['score'] = 1


# In[89]:


original_data2 = pd.pivot_table(original_data, index='orf', columns='Phenotype', values='score')


# In[90]:


original_data2[original_data2.isnull()] = 0


# In[91]:


original_data2.columns


# In[92]:


original_data2['Groups I, VI'] = original_data2[['Groups  I ,VI','Groups I, VI']].sum(axis=1)
original_data2.drop(columns=['Groups  I ,VI'], inplace=True)


# In[93]:


original_data2['data1'] = original_data2['Group I']*(-0.5) + original_data2['Group II'] * (0.5) + original_data2['Group III'] * (-1) + original_data2['Groups I, VI'] * (-0.5)
original_data2['data2'] = original_data2['Group III']*(-1) + original_data2['Group IV'] * (-0.5) + original_data2['Group V'] * (0.5)
original_data2['data3'] = original_data2['Group VI'] * 1 + original_data2['Groups I, VI'] * 1


# In[94]:


original_data2.head()


# In[95]:


original_data2 = original_data2[['data1','data2','data3']]


# In[96]:


original_data2 = original_data2.groupby(original_data2.index).mean()


# In[97]:


original_data2.shape


# In[98]:


original_data2[original_data2>1] = 1


# # Load & process tested strains

# In[99]:


tested = pd.read_excel('raw_data/ Gene List_ResGenLibrary.xlsx', sheet_name='mat_a_101501.txt', skiprows=1)


# In[100]:


tested.head()


# In[101]:


tested['orf'] = tested['ORF name'].astype(str)


# In[102]:


tested['orf'] = clean_orf(tested['orf'])


# In[103]:


tested.loc[tested['orf']=='YLR287-A','orf'] = 'YLR287C-A'


# In[104]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[105]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[106]:


tested_orfs = tested['orf'].unique()


# In[107]:


missing = [orf for orf in original_data2.index.values if orf not in tested_orfs]
missing


# In[108]:


original_data2 = original_data2.reindex(index=tested_orfs, fill_value=0)


# In[109]:


original_data2.head()


# # Prepare the final dataset

# In[110]:


data = original_data2.copy()


# In[111]:


dataset_ids = [11828, 11829, 5180]
datasets = datasets.reindex(index=dataset_ids)


# In[112]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[113]:


data.head()


# ## Subset to the genes currently in SGD

# In[114]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[115]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[116]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[117]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[118]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[119]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[120]:


from IO.save_data_to_db3 import *


# In[121]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




