#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 22615293
paper_name = 'hoepfner_parker_2012' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


compounds = ['Cpd1','Cpd2','Cpd3','Cpd4','Vori']


# In[15]:


original_data_list = {'HOP': [], 'HIP': []}
for c in compounds:
    original_data = pd.read_csv('raw_data/HIP-HOP Scores ' + c + '.txt', sep='\t')
    
    for et in ['HOP','HIP']:
        original_data1 = original_data.loc[original_data['EXPERIMENT_TYPE']==et,:].copy()
        original_data11 = pd.pivot_table(original_data1, index='GENE_NAME', columns='COMPOUND_CONCENTRATION', values='SCORE')
    
        original_data11 = original_data11.reset_index()
        original_data11['GENE_NAME'] = clean_genename(original_data11['GENE_NAME'])
        original_data11['orf'] = translate_sc(original_data11['GENE_NAME'], to='orf')
        
        # Make sure everything translated ok
        t = looks_like_orf(original_data11['orf'])
#         print(original_data11.loc[~t,])
        original_data11 = original_data11.loc[t,:]
        
        original_data11.set_index('orf', inplace=True)
        original_data11 = original_data11.groupby(original_data11.index).mean()
        
        # Rename columns
        cols = original_data11.columns.values
        cols = [('%s_%s' % (c, d.split('_')[1])) for d in cols]
        original_data11.columns = cols
        
        print(original_data11.shape)
        print(cols)
        
        original_data_list[et].append(original_data11)


# In[16]:


original_data1 = pd.concat(original_data_list['HOP'], axis=1)
original_data2 = pd.concat(original_data_list['HIP'], axis=1)


# In[17]:


doses = pd.read_excel('raw_data/doses_datasetids.xlsx', sheet_name='Sheet1', header=None)


# In[18]:


doses[4] = doses[2].apply(lambda x: ('%f' % x).rstrip('0').rstrip('.'))


# In[20]:


doses[5] = doses[0] + '_' + doses[4].astype(str)


# In[23]:


doses1 = doses.loc[doses[1]=='HOP',:].copy()
doses2 = doses.loc[doses[1]=='HIP',:].copy()


# In[28]:


doses1.set_index(5, inplace=True)
doses2.set_index(5, inplace=True)


# In[30]:


doses1 = doses1.reindex(index=original_data1.columns.values)
doses2 = doses2.reindex(index=original_data2.columns.values)


# In[34]:


original_data1.columns = doses1[3].values
original_data2.columns = doses2[3].values


# In[35]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[37]:


original_data.columns


# In[44]:


original_data.index.name = 'orf'


# # Prepare the final dataset

# In[45]:


data = original_data.copy()


# In[46]:


dataset_ids = original_data.columns.values
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




