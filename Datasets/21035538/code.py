#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 21035538
paper_name = 'uluisik_koc_2011' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[97]:


files = ['mmc3.xlsx','mmc2.xlsx']


# In[98]:


data_switch = {'+': 1,'-': 0}


# In[99]:


original_data_list = []
for f  in files:
    original_data = pd.read_excel('raw_data/'+f, sheet_name='Sayfa1', skiprows=2)
    print('Original data dimensions: %d x %d' % (original_data.shape))
    print(original_data.head())
    original_data['orf'] = original_data['Mutants'].astype(str)
    # Eliminate all white spaces & capitalize
    original_data['orf'] = clean_orf(original_data['orf'])
    # Translate to ORFs 
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    # Make sure everything translated ok
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    
    wt = original_data.iloc[0,:]
    
    original_data = original_data.loc[t,:]
    
    original_data.set_index('orf', inplace=True)
    original_data.drop(columns=['Mutants'], inplace=True)
    
    for c in original_data.columns:
        original_data[c] = original_data[c].apply(lambda x: data_switch[x])
    
    wt = wt[original_data.columns]
    wt = wt.apply(lambda x: data_switch[x])
    
    original_data = original_data.groupby(original_data.index).mean()
        
    # Subtract WT
    original_data = original_data.sub(wt, axis=1)
    
    # Subtract 0 nM
    original_data = original_data.sub(original_data['0 mM'], axis=0)
    
    original_data.drop(columns = ['0 mM'], inplace=True)
    
    print(original_data.shape)
    original_data_list.append(original_data)


# In[100]:


original_data = pd.concat(original_data_list, axis=1)


# In[101]:


original_data[original_data.isnull()] = 0


# In[102]:


original_data.head()


# In[103]:


original_data.index.name = 'orf'


# In[104]:


original_data.shape


# # Load & process tested strains

# In[105]:


tested = pd.read_csv('raw_data/tested_strains.txt', sep='\t', header=None)


# In[106]:


tested.head()


# In[107]:


tested['orf'] = tested[0].astype(str)


# In[108]:


tested['orf'] = clean_orf(tested['orf'])


# In[109]:


tested.loc[tested['orf']=='YLR287-A','orf'] = 'YLR287C-A'


# In[110]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[111]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[112]:


tested_orfs = tested['orf'].unique()


# In[113]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[114]:


tested_orfs = list(tested_orfs) + missing


# In[115]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[116]:


data = original_data.copy()


# In[117]:


dataset_ids = [438, 439, 440, 441, 442, 443, 147, 444, 445]
datasets = datasets.reindex(index=dataset_ids)


# In[118]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[119]:


data.head()


# ## Subset to the genes currently in SGD

# In[120]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[121]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[122]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[123]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[124]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[125]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[126]:


from IO.save_data_to_db3 import *


# In[127]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




