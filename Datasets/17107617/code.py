#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 17107617
paper_name = 'freimoser_amrhein_2006' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[7]:


mp = pd.read_excel('raw_data/Strains_plates_annotation.xlsx', sheet_name='Sheet1', skiprows=2)
mp.head()


# In[9]:


mp = mp.loc[mp['col'].notnull(),]


# In[10]:


mp['Plate'] = mp['Plate'].astype(int)
mp['col'] = mp['col'].astype(int)


# In[20]:


mp = mp[['Plate','row','col','ORF']].copy()
mp.columns = ['plate','row','col','orf']


# In[17]:


c = pd.read_csv('raw_data/columns_to_extract.txt', sep='\t')
c.head()


# In[18]:


# Transform to numbers:
alph = 'ABCDEFGHIJKLMNOPQRSTUVWXYZ'
c['Plate'] = c['Plate'].apply(lambda x: alph.index(x))
c['Row'] = c['Row'].apply(lambda x: alph.index(x))
c['Col'] = c['Col'].apply(lambda x: alph.index(x))
c['Data'] = c['Data'].apply(lambda x: alph.index(x))


# In[19]:


c.head()


# In[31]:


files = ['Plates1-10.xlsx','Plates11-20.xlsx','Plates21-30.xlsx','Plates31-40.xlsx','Plates41-50.xlsx','Plate51.xlsx']


# In[32]:


original_data.head()


# In[29]:


plate_col


# In[35]:


original_data_list = []
for f in files:
    xl = pd.ExcelFile('raw_data/' + f)
    
    sheet_names = c.loc[c['File']==f,'Sheet'].values
    
    for s in sheet_names:
        original_data = xl.parse(s, header=None)
        
        plate_col = c.loc[(c['File']==f) & (c['Sheet']==s),'Plate'].values[0]
        row_col = c.loc[(c['File']==f) & (c['Sheet']==s),'Row'].values[0]
        col_col = c.loc[(c['File']==f) & (c['Sheet']==s),'Col'].values[0]
        data_col = c.loc[(c['File']==f) & (c['Sheet']==s),'Data'].values[0]
        
        print('Original data dimensions: %d x %d' % (original_data.shape))
#         print(original_data.head())
        
        original_data = original_data.iloc[:, [plate_col, row_col, col_col, data_col]].copy()
        original_data.columns = ['plate','row','col','data']
        
        original_data = original_data.loc[original_data['plate'].notnull()]
        for cl in ['plate','col']:
            original_data[cl] = original_data[cl].astype(int)
        
        original_data = original_data.merge(mp, how='left', on=['plate','row','col'])
        
        original_data['orf'] = original_data['orf'].astype(str)
        original_data['orf'] = clean_orf(original_data['orf'])
        
        original_data['orf'] = translate_sc(original_data['orf'], to='orf')
        t = looks_like_orf(original_data['orf'])
        print(original_data.loc[~t,])
        
        original_data['data'] = pd.to_numeric(original_data['data'], errors='coerce')
        original_data.set_index('orf', inplace=True)
        original_data = original_data[['data']].copy()
        original_data = original_data.groupby(original_data.index).mean()
        
        print(original_data.shape)
        
        original_data_list.append(original_data)


# In[36]:


original_data = pd.concat(original_data_list, axis=0)


# In[37]:


original_data.shape


# In[38]:


original_data.head()


# In[44]:


t = looks_like_orf(original_data.index.values)
print(original_data.loc[~np.array(t),])


# In[45]:


original_data = original_data.loc[np.array(t),:]


# In[46]:


original_data = original_data.groupby(original_data.index).mean()


# In[47]:


original_data.shape


# # Prepare the final dataset

# In[48]:


data = original_data.copy()


# In[49]:


dataset_ids = [119]
datasets = datasets.reindex(index=dataset_ids)


# In[50]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[51]:


data.head()


# ## Subset to the genes currently in SGD

# In[52]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[53]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[54]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[55]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[56]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[57]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[58]:


from IO.save_data_to_db3 import *


# In[59]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




