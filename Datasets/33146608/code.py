#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 33146608
paper_name = 'kintaka_moriya_2020' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


sga_results_files = ['elife-54080-fig2-data1-v2.xlsx','Data_S3_tGFP_SGA_raw data.xlsx','Data_S4_NES-tGFP_SGA_raw data.xlsx']


# In[6]:


sheet_names = [['DMA U Rep1','DMA U Rep2'],['DMA LU Rep1','DMA LU Rep2']]
original_data_list = []


# In[7]:


dataset_ids = [16674,16675,16681,16677,16682,16678,16683,16685,16684]


# In[8]:


num_dataset = 0
for f in sga_results_files:
    print('File %s' % f)
    for sheet_names_this in sheet_names:
        this_data = pd.DataFrame()
        for s in sheet_names_this:
            print('Sheet %s' % s)
            t = pd.read_excel('raw_data/' + f, sheet_name=s)
            print('Original data dimensions: %d x %d' % (t.shape))
            t = t[['Array ORF','Score']]
            t['Score'] = t['Score'].astype(float)
            this_data = pd.concat((this_data, t), axis=0)
        this_data.columns = ['orf', dataset_ids[num_dataset]]
        num_dataset += 1
        original_data_list.append(this_data)


# In[9]:


gfp_expression_file = 'elife-54080-fig4-data1-v2.xlsx'
sheet_names = ['GFPunit (GFP)','GFPunit(NES-tGFP)','GFPunit(tGFP)']


# In[10]:


for s in sheet_names:
    t = pd.read_excel('raw_data/' + gfp_expression_file, sheet_name=s)
    print('Original data dimensions: %d x %d' % (t.shape))
    t = t.loc[t['array name']=='DMA']
    t = t[['ORF','GFPunit_Average']]
    t['GFPunit_Average'] = t['GFPunit_Average'].astype(float)
    t.columns = ['orf', dataset_ids[num_dataset]]
    num_dataset += 1
    original_data_list.append(t)


# In[11]:


len(original_data_list)


# In[12]:


data_list = []
for df in original_data_list:
    df['orf'] = df['orf'].astype(str)
    df['orf'] = clean_orf(df['orf'])
    df = df.groupby('orf').mean().reset_index()
    df['orf'] = translate_sc(df['orf'], to='orf')
    t = looks_like_orf(df['orf'])
    print(df.loc[~t,])
    df = df.groupby('orf').mean()
    data_list.append(df)


# In[13]:


data = pd.concat(data_list, axis=1, join='outer')


# In[14]:


data.shape


# In[15]:


data.index.name='orf'


# # Prepare the final dataset

# In[16]:


datasets = datasets.reindex(index=dataset_ids)


# In[17]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[18]:


data.head()


# ## Subset to the genes currently in SGD

# In[19]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[20]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])


# In[21]:


data.head()


# # Normalize

# In[22]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[23]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[24]:


data_norm[data.isnull()] = np.nan


# In[25]:


data_all = data.join(data_norm)


# In[26]:


data_all.head()


# # Print out

# In[27]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[28]:


from IO.save_data_to_db3 import *


# In[29]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




