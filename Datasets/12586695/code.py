#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 12586695
paper_name = 'enyenihi_saunders_2003' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[58]:


original_data = pd.read_excel('raw_data/Table_1_of_total_results_repaired.xlsx', sheet_name='Sheet1', header=None)


# In[59]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[60]:


original_data.head()


# In[61]:


original_data_list = []
for c in np.arange(0,30,2):
    t = original_data[[c,c+1]].copy()
    t.columns = [0,1]
    original_data_list.append(t)


# In[62]:


original_data = pd.concat(original_data_list, axis=0, ignore_index=True)


# In[63]:


original_data.head()


# In[64]:


original_data = original_data.loc[original_data[0].notnull()]


# In[65]:


original_data['gene'] = original_data[0].astype(str)


# In[66]:


# Eliminate all white spaces & capitalize
original_data['gene'] = clean_genename(original_data['gene'])


# In[67]:


typo_fixes = {'YBLO29C':'YBL029C','HDR1':'YBR138C','PET1OO':'YDR079W','SRI1':'YEL025C','TOS9':'YEL007W','TOS10':'YGR153W',
              'KRE20': 'YAL056C-A','KRE23':'YAL042C-A','FYV14':'YDL213C','FYV9':'YDR140W','FYV11':'YFL023W','FYV2':'YIL054W','EFR4':'YLR114C',
              'RPSOA':'YGR214W','GIF1':'YIR024C','HSP-150':'YJL159W','GLGL1':'YKR058W','TOS5':'YKR011C','TOS7':'YOL019W','WH12':'YOR043W'}
original_data['gene'] = original_data['gene'].apply(lambda x: typo_fixes[x] if x in typo_fixes.keys() else x)


# In[68]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['gene'].values, to='orf')


# In[71]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[72]:


original_data = original_data.loc[t,:]


# In[74]:


# A little workaround
original_data['data'] = clean_genename(original_data[1])


# In[77]:


data_switch= {'NORMAL':0, 'L0W':-1, 'LOW': -1, 'LOW-4': -1, 'LOWNO-TETRADS':-1, 'LOWNOTETRAD':-1, 'LOW-DYADS': -1, 'VLOW': -2, 'VLOW-4': -2, 'NONE': -3}
original_data['data'] = original_data['data'].apply(lambda x: data_switch[x])


# In[78]:


original_data.set_index('orf', inplace=True)


# In[79]:


original_data = original_data[['data']].copy()


# In[80]:


original_data = original_data.groupby(original_data.index).mean()


# In[81]:


original_data.shape


# # Prepare the final dataset

# In[82]:


data = original_data.copy()


# In[83]:


dataset_ids = [72]
datasets = datasets.reindex(index=dataset_ids)


# In[84]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[85]:


data.head()


# ## Subset to the genes currently in SGD

# In[86]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[87]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[88]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[89]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[90]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[91]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[92]:


from IO.save_data_to_db3 import *


# In[93]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




