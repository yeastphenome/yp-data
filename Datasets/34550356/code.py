#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 34550356
paper_name = 'zhou_foijer_2021' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[111]:


sheets = ['Screen1-YKO raw data','Screen2-YKO raw data']


# In[112]:


original_data = []
for s in sheets:
    original_data1 = pd.read_excel('raw_data/FileS1.xls', sheet_name=s, skiprows=16)
    print('Original data dimensions: %d x %d' % (original_data1.shape))
    original_data1['orf'] = original_data1['ID Column'].astype(str)
    # Eliminate all white spaces & capitalize
    original_data1['orf'] = clean_orf(original_data1['orf'])
    # Translate to ORFs 
    original_data1['orf'] = translate_sc(original_data1['orf'], to='orf')
    # Make sure everything translated ok
    t = looks_like_orf(original_data1['orf'])
    print(original_data1.loc[~t,])
    original_data1 = original_data1.loc[t,]
    
    t = original_data1['Normalized Growth Ratio (Comparer::Exp)'].str.contains('excluded')
    original_data1 = original_data1.loc[~t,]
    
    original_data1['data'] = 1/pd.to_numeric(original_data1['Growth Ratio (Comparer / Exp)'], errors='coerce')
    original_data1 = pd.pivot_table(original_data1, index='orf', columns='Query', values='data')
    original_data1 = original_data1.groupby(original_data1.index).mean()
    print(original_data1.shape)
    
    original_data.append(original_data1)


# In[113]:


sheets = ['Screen3-MELK_YKO raw data','Screen3-MELK_KD_YKO raw data']
data_cols = ['MELK_YKO3_size.mean.norm_div_control','MELK_KD_YKO3_size.mean.norm_div_control']
for ix_s, s in enumerate(sheets):
    original_data2 = pd.read_excel('raw_data/FileS1.xls', sheet_name=s)
    print('Original data dimensions: %d x %d' % (original_data2.shape))
    original_data2['orf'] = original_data2['ORF'].astype(str)
    # Eliminate all white spaces & capitalize
    original_data2['orf'] = clean_orf(original_data2['orf'])
    # Translate to ORFs 
    original_data2['orf'] = translate_sc(original_data2['orf'], to='orf')
    # Make sure everything translated ok
    t = looks_like_orf(original_data2['orf'])
    print(original_data2.loc[~t,])
    
    original_data2 = original_data2.loc[t,]
    original_data2['data'] = pd.to_numeric(original_data2[data_cols[ix_s]], errors='coerce')
    original_data2.set_index('orf', inplace=True)
    original_data2 = original_data2[['data']]
    original_data2 = original_data2.groupby(original_data2.index).mean()
    print(original_data2.shape)
    
    original_data2.columns = [s]
    
    original_data.append(original_data2)


# In[114]:


len(original_data)


# In[162]:


original_data_all = pd.concat(original_data, axis=1)


# In[163]:


original_data_all.head()


# In[164]:


original_data_all.columns = ['MELK','MELK_kd','MELK_kd','MELK','MELK','MELK_kd']


# In[165]:


original_data_all = original_data_all.T
original_data_all = original_data_all.groupby(original_data_all.index.values).mean()
original_data_all = original_data_all.T


# In[166]:


original_data_all.head()


# In[167]:


original_data_all.shape


# # Prepare the final dataset

# In[168]:


data = original_data_all.copy()


# In[169]:


data.index.names = ['orf']


# In[170]:


dataset_ids = [22059,22060]
datasets = datasets.reindex(index=dataset_ids)


# In[171]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[172]:


data.head()


# ## Subset to the genes currently in SGD

# In[173]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[174]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[175]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[176]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[177]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# In[178]:


data_all.loc[(slice(None),'YOL025W'),:]


# # Print out

# In[179]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[180]:


from IO.save_data_to_db3 import *


# In[181]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




