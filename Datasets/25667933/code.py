#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 25667933
paper_name = 'nislow_hammond_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[6]:


files = ['TableS5_vsT1-hom1-dropsToBg.xlsx',
         'TableS6_vsT1-hom1-NaCl-dropsToBg.xlsx',
         'TableS8_vsT1-het1-dropsToBg.xlsx',
         'TableS9_vsT1-het1-NaCl-dropsToBg.xlsx']
sheet_names = ['ground','flight']


# In[8]:


original_data_list = []
for f in files:
    for s in sheet_names:
        original_data = pd.read_excel('raw_data/' + f, sheet_name=s)
        
        original_data['orf'] = original_data['Strain'].apply(lambda x: x.split(':')[0])
        original_data['orf'] = clean_orf(original_data['orf'])
        original_data['orf'] = translate_sc(original_data['orf'], to='orf')
        # Make sure everything translated ok
        t = looks_like_orf(original_data['orf'])
        print(original_data.loc[~t,])
        original_data.set_index('orf', inplace=True)
        
        data_cols = [x for x in original_data.columns.values if 'Log2ratio' in x]
        original_data = original_data.loc[:,data_cols].copy()
        original_data = original_data.groupby(original_data.index).mean()
        
        print(original_data.shape)
        original_data_list.append(original_data)


# In[32]:


original_data = pd.concat(original_data_list, axis=1)


# In[33]:


original_data.index.name='orf'


# In[34]:


# Data represent log(gen1/gen2) such that, originally, higher values correspond to more extreme depletion
original_data = -original_data


# In[35]:


original_data.head()


# # Prepare the final dataset

# In[36]:


data = original_data.copy()


# In[37]:


dataset_ids = [5276, 5277, 5278, 5279, 5280, 5281, 5282, 5283, 5284, 5285, 5286, 5287, 5288, 5289]
datasets = datasets.reindex(index=dataset_ids)


# In[38]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[39]:


data.head()


# ## Subset to the genes currently in SGD

# In[40]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[41]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[42]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[43]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[44]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[45]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[46]:


from IO.save_data_to_db3 import *


# In[47]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




