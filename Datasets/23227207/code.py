#!/usr/bin/env python
# coding: utf-8

# In[4]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[5]:


paper_pmid = 23227207
paper_name = 'serviene_urbonavicius_2012' 


# In[6]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[7]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[60]:


files = ['pone.0050779.s004.xls','pone.0050779.s005.xls']


# In[61]:


def get_score(s):
    s = str(s).replace(" ","")
    scale = {'+/-':1,'+':2,'++':3,'+++':4}
    
    if s[0] == 'R':
        score = scale[s[1:]]
    elif s[0] == 'S':
        score = -scale[s[1:]]
    else:
        score = np.nan
    return score


# In[62]:


original_data_list = []

for f  in files:
    original_data = pd.read_excel('raw_data/' + f, sheet_name='Sheet1', skiprows=10)
    print('Original data dimensions: %d x %d' % (original_data.shape))
    original_data['orf'] = original_data.iloc[:,0].astype(str)
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    t = looks_like_orf(original_data['orf'])
    print(original_data.loc[~t,])
    original_data = original_data.loc[t,:]
    
    for f in ['Unnamed: 3','Unnamed: 4','Unnamed: 5']:
        original_data[f] = original_data[f].apply(get_score)
        
    original_data.set_index('orf', inplace=True)
    original_data = original_data.loc[:,['Unnamed: 3','Unnamed: 4','Unnamed: 5']]
    original_data = original_data.groupby(original_data.index).mean()
    
    print(original_data.shape)
    original_data_list.append(original_data)


# In[63]:


original_data1, original_data2 = original_data_list


# In[64]:


original_data = original_data1.join(original_data2, how='outer', lsuffix='_1', rsuffix='_2')


# In[66]:


original_data['data'] = original_data.sum(axis=1)


# In[69]:


original_data.sort_values(by='data', ascending=False).head()


# # Prepare the final dataset

# In[70]:


data = original_data[['data']].copy()


# In[71]:


dataset_ids = [16525]
datasets = datasets.reindex(index=dataset_ids)


# In[72]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[73]:


data.head()


# ## Subset to the genes currently in SGD

# In[74]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[75]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[76]:


data.head()


# # Normalize

# In[77]:


data_norm = normalize_phenotypic_scores(data, has_tested=False)


# In[78]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[79]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)


# In[80]:


data_all.head()


# # Print out

# In[81]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[82]:


from IO.save_data_to_db3 import *


# In[83]:


save_data_to_db(data_all, paper_pmid)


# In[84]:


data_all.shape


# In[ ]:




