#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 26341223
paper_name = 'garcia_arroyo_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['pmid', 'name'])


# In[4]:


datasets.set_index('pmid', inplace=True)


# # Load & process the data

# ### Zymolyase

# In[5]:


path = 'raw_data/screening zymo/'
excel_files = [f for f in os.listdir(path) if os.path.isfile(os.path.join(path, f))]


# In[6]:


files_to_sheets = pd.read_excel('raw_data/zymo_files_to_sheets.xlsx', sheet_name='Sheet1', header=None)


# In[7]:


all_data = pd.DataFrame()
for f in excel_files:
    sheet_name = files_to_sheets.loc[files_to_sheets[0]==f,1].values[0]
    d = pd.read_excel(os.path.join(path, f), sheet_name=sheet_name)
    d.columns = [c.lower() for c in d.columns]
    all_data = pd.concat([all_data, d], axis=0)


# In[8]:


all_data.drop(columns=['unnamed: 3'], inplace=True)


# In[9]:


all_data['orf'] = all_data['orf'].astype(str)


# In[10]:


# Eliminate all white spaces & capitalize
all_data['orf'] = clean_orf(all_data['orf'])


# In[11]:


all_data = all_data.groupby('orf').mean()


# In[12]:


all_data = all_data.reset_index()


# In[13]:


# Translate to ORFs 
all_data['orf'] = translate_sc(all_data['orf'], to='orf')


# In[14]:


# Make sure everything translated ok
t = looks_like_orf(all_data['orf'])
print(all_data.loc[~t,])


# In[15]:


all_data.set_index('orf', inplace=True)


# In[16]:


all_data = all_data.div(all_data.loc['WT',:])


# In[17]:


all_data.drop(index='WT', inplace=True)


# In[18]:


all_data = all_data['24h'].to_frame()
all_data.columns = ['ZYM sens']


# In[19]:


all_data = all_data.groupby(all_data.index).mean()


# In[20]:


all_data.shape


# ### Congo Red and caspofungin

# In[21]:


all_data2 = pd.read_excel('raw_data/12864_2015_1879_MOESM2_ESM-2.xls', sheet_name='Hoja1', skiprows=1)


# In[22]:


all_data2['CR_num'] = all_data2['CR'].apply(lambda x: len(x) if isinstance(x, str) else x)


# In[23]:


all_data2['CAS_num'] = all_data2['CAS'].apply(lambda x: len(x) if isinstance(x, str) else x)


# In[24]:


all_data2['ZYM_num'] = all_data2['ZYM'].apply(lambda x: len(x) if isinstance(x, str) else x)


# In[25]:


all_data2 = all_data2[['ORF','CR_num','CAS_num','ZYM_num']]


# In[26]:


all_data2['ORF'] = all_data2['ORF'].astype(str)


# In[27]:


# Eliminate all white spaces & capitalize
all_data2['ORF'] = clean_orf(all_data2['ORF'])


# In[28]:


# Translate to ORFs 
all_data2['ORF'] = translate_sc(all_data2['ORF'], to='orf')


# In[29]:


# Make sure everything translated ok
t = looks_like_orf(all_data2['ORF'])
print(all_data2.loc[~t,])


# In[30]:


all_data2 = all_data2.loc[t,]


# In[31]:


all_data2 = all_data2.groupby('ORF').mean()


# In[32]:


all_data2 = -all_data2


# ### Caspofungin resistance

# In[33]:


all_data3 = pd.read_excel('raw_data/Table3.xlsx', sheet_name='Sheet1')


# In[34]:


all_data3['ORF'] = all_data3['ORF'].astype(str)


# In[35]:


# Eliminate all white spaces & capitalize
all_data3['ORF'] = clean_orf(all_data3['ORF'])


# In[36]:


# Translate to ORFs 
all_data3['ORF'] = translate_sc(all_data3['ORF'], to='orf')


# In[37]:


# Make sure everything translated ok
t = looks_like_orf(all_data3['ORF'])
print(all_data3.loc[~t,])


# In[38]:


all_data3 = all_data3[['ORF']].copy()


# In[39]:


all_data3['CAS res'] = 1


# In[40]:


all_data3 = all_data3.groupby('ORF').mean()


# # Merge all datasets

# In[64]:


original_data = all_data.join(all_data2, how='outer')
original_data = original_data.join(all_data3, how='outer')


# In[65]:


# Dropping the discrete phenotypes for zymolyase because it is replaced by the quantitative data in ZYM sens
original_data.drop(columns=['ZYM_num'], inplace=True)


# In[66]:


# Set all NaN values to 0, effectively assuming that the quantitative ZYM sens dataset represents the tested universe for the other experiments as well
for c in [1,2,3]:
    col = original_data.columns.values[c]
    original_data.loc[original_data.loc[:,col].isnull(), col] = 0


# In[67]:


original_data.index.name = 'orf'


# # Prepare the final dataset

# In[68]:


data = original_data.copy()


# In[69]:


dataset_ids = [16465,16466,16467,16470]
datasets = datasets.reindex(index=dataset_ids)


# In[70]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[71]:


data.head()


# ## Subset to the genes currently in SGD

# In[72]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[73]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[74]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[75]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[76]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[77]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[78]:


from IO.save_data_to_db3 import *


# In[79]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




