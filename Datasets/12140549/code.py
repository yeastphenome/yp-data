#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 12140549
paper_name = 'giaever_johnston_2002' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# In[5]:


datasets.shape


# # Load & process the data

# In[6]:


data_files = [f for f in os.listdir('raw_data/') if f.endswith('.txt') and not f.startswith('ypd') and not f.startswith('Cell')]
len(data_files)


# In[7]:


original_data_list = []
original_data_experiment_list = []
for ixf, f in enumerate(data_files):
    print(ixf)
    experiment_id = int(f.split('_')[0])
    sign = -1 if f.split('_')[-1] == 'sen.txt' else 1
    original_data = pd.read_csv('raw_data/' + f, header=None, sep='\t')
#     print('Original data dimensions: %d x %d' % (original_data.shape))
#     print(original_data.head())
    original_data['orf'] = original_data[0].astype(str)    
    original_data['orf'] = clean_orf(original_data['orf'])
    original_data['orf'] = translate_sc(original_data['orf'], to='orf')
    t = looks_like_orf(original_data['orf'])
#     print(original_data.loc[~t,])
    
    original_data.set_index('orf', inplace=True)
    original_data['data'] = sign * pd.to_numeric(original_data[1], errors='coerce')
    original_data = original_data[['data']].copy()
    original_data = original_data.groupby(original_data.index).mean()
    
    original_data_list.append(original_data)
    original_data_experiment_list.append(experiment_id)


# In[8]:


original_data = pd.concat(original_data_list, axis=1)


# In[9]:


original_data.columns = original_data_experiment_list


# In[10]:


original_data.shape


# In[11]:


# Load dataset ids
dt = pd.read_excel('raw_data/phenotype_mapping.xlsx', sheet_name='Sheet1')


# In[12]:


dt.set_index('Experiment', inplace=True)


# In[13]:


dt = dt.reindex(index=original_data_experiment_list)


# In[14]:


dataset_ids = dt['Dataset id'].values


# In[15]:


dataset_ids


# In[16]:


original_data.columns = dataset_ids


# In[17]:


original_data = original_data.T
original_data = original_data.groupby(original_data.index).mean()
original_data = original_data.T


# In[18]:


original_data.shape


# In[19]:


original_data.notnull().sum(axis=0)


# # Load data (2)

# In[20]:


ypd = pd.read_csv('raw_data/ypd.txt', sep='\t')


# In[21]:


ypd.head()


# In[22]:


ypd['orf'] = ypd['ORF'].astype(str)    
ypd['orf'] = clean_orf(ypd['orf'])
ypd['orf'] = translate_sc(ypd['orf'], to='orf')
t = looks_like_orf(ypd['orf'])
print(ypd.loc[~t,])


# In[23]:


ypd.set_index('orf', inplace=True)
ypd['data'] = -pd.to_numeric(ypd['Average Ratio'], errors='coerce')
ypd = ypd[['data']].copy()
ypd = ypd.groupby(ypd.index).mean()


# In[24]:


ypd.head()


# In[25]:


ypd.columns = [16187]


# In[26]:


ypd.shape


# In[27]:


original_data2 = original_data.join(ypd, how='outer')


# In[28]:


original_data2.shape


# In[30]:


# Set missing YPD values to 0
original_data2.loc[original_data2[16187].isnull(), 16187] = 0


# # Load data (3)

# In[31]:


morph = pd.read_csv('raw_data/Cell_Morph_Screen_Table.txt', sep='\t')


# In[32]:


morph.head()


# In[33]:


morph.columns = [c.strip() for c in morph.columns]


# In[34]:


morph['orf'] = morph['ORF'].astype(str)    
morph['orf'] = clean_orf(morph['orf'])
morph.loc[morph['orf']=='YELOO1C','orf'] = 'YEL001C'
morph['orf'] = translate_sc(morph['orf'], to='orf')
t = looks_like_orf(morph['orf'])
print(morph.loc[~t,])
morph = morph.loc[t,:]


# In[35]:


mp = pd.read_excel('raw_data/phenotype_mapping2.xlsx', sheet_name='Sheet1')
mp.head()


# In[36]:


mp.set_index('ORIGINAL', inplace=True)


# In[37]:


for d in mp['Dataset id'].unique():
    morph[d] = 0


# In[38]:


for ixr, row in morph.iterrows():
    ps = [x.strip() for x in row['Cell Shape Morphologies'].split(';')]
    
    for p in ps:
        parts = p.split(' ')
        if len(parts) > 1:
            ph = parts[0]
            try:
                score = int(parts[1])
            except ValueError as e:
                next
    
            if ph in mp.index.values:
                morph.loc[ixr, mp.loc[ph,'Dataset id']] = score * mp.loc[ph, 'COEFFICIENT']


# In[39]:


morph.set_index('orf', inplace=True)
morph = morph[[725, 726, 727, 729, 728]].copy()


# In[40]:


morph = morph.groupby(morph.index).mean()


# In[41]:


morph.sum(axis=0)


# In[42]:


original_data2 = original_data2.join(morph, how='outer')


# In[43]:


original_data2.shape


# In[44]:


original_data2.index.name='orf'


# In[45]:


original_data2.columns


# In[46]:


original_data2.notnull().sum(axis=0)


# # Prepare the final dataset

# In[47]:


data = original_data2.copy()


# In[48]:


dataset_ids = original_data2.columns.values
datasets = datasets.reindex(index=dataset_ids)


# In[49]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[50]:


data.head()


# ## Subset to the genes currently in SGD

# In[51]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[52]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[53]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[54]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[55]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[56]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[57]:


from IO.save_data_to_db3 import *


# In[58]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




