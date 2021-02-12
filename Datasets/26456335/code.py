#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')

import itertools


# # Initial setup

# In[2]:


paper_pmid = 26456335
paper_name = 'mccormick_kennedy_2015' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


original_data = pd.read_excel('raw_data/rls-summary-for-Anastasia-Baryshnikova-all-BY-haploid-deletion-YPD-30C-mm042018.xlsx', 
                            sheet_name='rls')


# In[6]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[7]:


# Fix a typo
original_data.loc[original_data['set_background']=='BY4,742','set_background'] = 'BY4742'
original_data['set_mating_type'] = original_data['set_mating_type'].str.lower()


# In[8]:


# Only keep BY4742 (systematic screen)
original_data = original_data.loc[original_data['set_background'] == 'BY4742',:]


# In[9]:


# Only keep single mutants
original_data['set_genotype'] = original_data['set_genotype'].str.strip()
original_data = original_data.loc[~original_data['set_genotype'].str.contains(' '),:]


# In[10]:


original_data.head()


# In[11]:


original_data['genes'] = original_data['set_genotype'].astype(str)


# In[12]:


# Eliminate all white spaces & capitalize
original_data['genes'] = clean_genename(original_data['genes'])


# In[13]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['genes'], to='orf')


# In[14]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,['set_name','genes','orf']])


# In[15]:


original_data.loc[~t,'orf'] = original_data.loc[~t,'set_name']


# In[16]:


manual_fixes = {'rpl20b':'YOR312C','sus1':'YBR111W-A','afg3::KanMX':'YER017C',
                'tor1':'YJR066W','pph22':'YDL188C','rpn4':'YDL020C',
                'scp1':'YOR367W','por1':'YNL055C','pmt3':'YOR321W','sir2':'YDL042C',
                'dbp3':'YGL078C','ymr226c': 'YMR226C'}

for typo in manual_fixes.keys():
    original_data.loc[original_data['orf']==typo,'orf'] = manual_fixes[typo]


# In[17]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,['set_name','genes','orf']])


# In[18]:


original_data = original_data.loc[t,:]


# In[19]:


original_data['rls'] = original_data['set_lifespans'].apply(lambda x: [int(t) for t in str(x).split(',')])


# In[20]:


original_data['ref_rls'] = original_data['ref_lifespans'].apply(lambda x: [int(t) for t in str(x).split(',')])


# In[21]:


# original_data.head()


# In[22]:


# Merge all raw measurements for all replicates
all_orfs = np.unique(original_data['orf'].values)
original_data2 = pd.DataFrame(index=all_orfs, columns=['rls','ref_rls'])


# In[23]:


all_orfs.shape


# In[24]:


for orf in all_orfs:
    this = original_data.loc[original_data['orf']==orf]
    original_data2.loc[orf,'rls'] = list(itertools.chain.from_iterable(this['rls']))
    original_data2.loc[orf,'ref_rls'] = list(itertools.chain.from_iterable(this['ref_rls']))


# In[25]:


original_data2['rls_num'] = original_data2['rls'].apply(lambda x: len(x))
original_data2['ref_rls_num'] = original_data2['ref_rls'].apply(lambda x: len(x))


# In[26]:


original_data2['rls_mean'] = original_data2['rls'].apply(lambda x: np.nanmean(np.array(x)))
original_data2['ref_rls_mean'] = original_data2['ref_rls'].apply(lambda x: np.nanmean(np.array(x)))


# In[27]:


original_data2['rls_ratio'] = original_data2['rls_mean'] / original_data2['ref_rls_mean']


# In[28]:


# % Only keep data with n > 5 (rest is unreliable)
# % The authors subsequently retested these strains but the raw version of that data is
# % (unfortunately) not recoverable. Only the published results.
# % Since all but one of the published (most reliable) strains has n > 5,
# % we've decided to set the n < 5 strains to 1 (instead of NaN) to indicate that they are
# % likely neither short-lived nor long-lived (instead of "not tested")

original_data2.loc[(original_data2['rls_num']<=5) & (original_data2['ref_rls_num']<=5),'rls_ratio'] = 1


# In[29]:


original_data2['data'] = original_data2['rls_ratio']


# In[30]:


original_data2 = original_data2[['data']].copy()


# In[31]:


original_data2.index.name='orf'


# In[32]:


original_data2.shape


# # Prepare the final dataset

# In[34]:


data = original_data2.copy()


# In[35]:


dataset_ids = [696]
datasets = datasets.reindex(index=dataset_ids)


# In[36]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[37]:


data.head()


# ## Subset to the genes currently in SGD

# In[38]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[39]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[40]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[41]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[42]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[43]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[44]:


from IO.save_data_to_db3 import *


# In[45]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




