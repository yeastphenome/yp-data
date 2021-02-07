#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 12096123
paper_name = 'wilson_roach_2002' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[8]:


original_data = pd.read_excel('raw_data/suppdata.xlsx', sheet_name='SuppDataREV', skiprows=1)


# In[9]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[10]:


original_data.head()


# In[11]:


original_data['orf'] = original_data['Unnamed: 2'].astype(str)


# In[12]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[16]:


original_data.loc[original_data['orf']=='YORO36W','orf'] = 'YOR036W'


# In[17]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[18]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[20]:


original_data = original_data.loc[t,:]


# In[22]:


data_switch = {'High': 1, 'High (pink)': 1, 'High (very)': 2, 'Low': -1, 'Low (very)': -2}
original_data['data'] = original_data['Unnamed: 1'].apply(lambda x: data_switch[x.strip()])


# In[23]:


original_data.set_index('orf', inplace=True)


# In[24]:


original_data = original_data[['data']].copy()


# In[25]:


original_data = original_data.groupby(original_data.index).mean()


# In[26]:


original_data.shape


# # Load & process tested strains

# In[29]:


tested = pd.read_excel('raw_data/ResGen Diploids inventory.xlsx', sheet_name='Inventory', skiprows=1)


# In[30]:


tested.head()


# In[31]:


tested['orf'] = tested['ORF name'].astype(str)


# In[32]:


tested['orf'] = clean_orf(tested['orf'])


# In[35]:


typos = {'TAL004W':'YAL004W','YELOO1C':'YEL001C','KL187C':'YKL187C'}
tested['orf'] = tested['orf'].apply(lambda x: typos[x] if x in typos.keys() else x)


# In[36]:


tested['orf'] = translate_sc(tested['orf'], to='orf')


# In[37]:


# Make sure everything translated ok
t = looks_like_orf(tested['orf'])
print(tested.loc[~t,])


# In[38]:


tested = tested.loc[t,:]


# In[39]:


tested_orfs = tested['orf'].unique()


# In[40]:


missing = [orf for orf in original_data.index.values if orf not in tested_orfs]
missing


# In[41]:


original_data = original_data.reindex(index=tested_orfs, fill_value=0)


# # Prepare the final dataset

# In[42]:


data = original_data.copy()


# In[43]:


dataset_ids = [4949]
datasets = datasets.reindex(index=dataset_ids)


# In[44]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[45]:


data.head()


# ## Subset to the genes currently in SGD

# In[46]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[47]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[48]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[49]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[50]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[51]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[52]:


from IO.save_data_to_db3 import *


# In[53]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




