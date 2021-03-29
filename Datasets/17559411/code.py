#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 17559411
paper_name = 'shuster_rosenberg_2007' 


# In[4]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[5]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[22]:


original_data = pd.read_excel('raw_data/femsyr_268_Table_S1.xlsx', sheet_name='Table 1', skiprows=1)


# In[23]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[24]:


original_data.head()


# In[25]:


for c in ['MACH negative','Reduced','Slow growers']:
    
    original_data[c] = original_data[c].apply(lambda x: [y.split('(')[0].strip().strip('\n') for y in str(x).split(',')])


# In[26]:


orf1 = [orf for list_ in original_data['MACH negative'] for orf in list_]


# In[28]:


orf2 = [orf for list_ in original_data['Reduced'] for orf in list_]


# In[30]:


orf3 = [orf for list_ in original_data['Slow growers'] for orf in list_]


# In[31]:


orfs = list(set(orf1+orf2+orf3))


# In[37]:


original_data = pd.DataFrame(index=orfs, columns=['mach','growth'], data=0)


# In[39]:


original_data.loc[orf1,'mach'] = -2
original_data.loc[orf1,'growth'] = 1
original_data.loc[orf2,'mach'] = -1
original_data.loc[orf2,'growth'] = 1
original_data.loc[orf3,'growth'] = -1


# In[41]:


original_data = original_data.reset_index()


# In[42]:


original_data.head()


# In[44]:


original_data['genename'] = original_data['index'].astype(str)


# In[45]:


# Eliminate all white spaces & capitalize
original_data['genename'] = clean_genename(original_data['genename'])


# In[46]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['genename'], to='orf')


# In[47]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[48]:


original_data = original_data.loc[t,:]


# In[49]:


original_data.set_index('orf', inplace=True)


# In[50]:


original_data = original_data[['mach','growth']].copy()


# In[51]:


original_data = original_data.groupby(original_data.index).mean()


# In[52]:


original_data.shape


# In[53]:


original_data.head()


# # Prepare the final dataset

# In[55]:


data = original_data.copy()


# In[56]:


dataset_ids = [16705, 21881]
datasets = datasets.reindex(index=dataset_ids)


# In[57]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[58]:


data.head()


# ## Subset to the genes currently in SGD

# In[59]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[60]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[61]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[62]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[63]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[64]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[65]:


from IO.save_data_to_db3 import *


# In[66]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




