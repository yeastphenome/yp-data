#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 15695358
paper_name = 'huang_oshea_2005' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load map

# In[46]:


original_data = pd.read_excel('raw_data/MATa Collection.xlsx', sheet_name='MATa Collection.xls', skiprows=3)


# In[47]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[48]:


original_data.head()


# In[49]:


original_data['orf'] = original_data['ORF name'].astype(str)


# In[50]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[51]:


original_data.loc[original_data['orf']=='YLR287-A','orf'] = 'YLR287C-A'


# In[52]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[53]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[54]:


original_data = original_data.loc[t,:]


# In[56]:


original_data = original_data.loc[original_data['column'].notnull() & original_data['row'].notnull() & original_data['plate'].notnull()]


# In[58]:


original_data['plate'] = original_data['plate'].astype(int)
original_data['column'] = original_data['column'].astype(int)


# In[61]:


original_data.set_index(['plate','row','column'], inplace=True)


# In[63]:


original_data = original_data[['orf']].copy()


# In[64]:


original_data.head()


# In[ ]:


original_data = original_data.groupby(original_data.index).mean()


# In[ ]:


original_data.shape


# # Load data

# In[65]:


dt = pd.read_excel('raw_data/Initial Screen Data.xls', sheet_name='Sheet1', skiprows=2)


# In[66]:


dt.head()


# In[67]:


dt['Unnamed: 0'] = dt['Unnamed: 0'].astype(str)


# In[68]:


ix_plates = dt.loc[dt['Unnamed: 0'].str.startswith('Plate')].index.values


# In[81]:


dt_this_list = []
for ix in ix_plates:
    plate_num = int(dt.loc[ix,'Unnamed: 0'].split('#')[1])
    dt_this = dt.loc[ix+1:ix+5,:].T
    
    dt_this.columns = ['pos','0','120','240','360']
    dt_this.drop(index='Unnamed: 0', inplace=True)
    
    dt_this['plate'] = plate_num
    dt_this['row'] = dt_this['pos'].apply(lambda x: x[0])
    dt_this['column'] = dt_this['pos'].apply(lambda x: int(x[1:]))

    dt_this.set_index(['plate','row','column'], inplace=True)
    dt_this.drop(columns=['pos'], inplace=True)
    
    dt_this_list.append(dt_this)


# In[82]:


dt = pd.concat(dt_this_list, axis=0)


# In[83]:


dt.shape


# In[87]:


original_data = dt.join(original_data, how='left')


# In[89]:


original_data = original_data.loc[original_data['orf'].notnull()]


# In[90]:


original_data.shape


# In[92]:


original_data.head()


# In[93]:


original_data.set_index('orf', inplace=True)


# In[95]:


original_data['data'] = original_data[['0','120','240','360']].mean(axis=1)


# In[96]:


original_data = original_data[['data']].copy()


# In[98]:


original_data = original_data.groupby(original_data.index).mean()


# In[99]:


original_data.head()


# In[100]:


original_data.shape


# # Prepare the final dataset

# In[101]:


data = original_data.copy()


# In[102]:


dataset_ids = [139]
datasets = datasets.reindex(index=dataset_ids)


# In[103]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[104]:


data.head()


# ## Subset to the genes currently in SGD

# In[105]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[106]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[107]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[108]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[109]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[110]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[111]:


from IO.save_data_to_db3 import *


# In[112]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




