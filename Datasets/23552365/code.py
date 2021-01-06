#!/usr/bin/env python
# coding: utf-8

# In[1]:


get_ipython().run_line_magic('run', '../yp_utils.py')


# # Initial setup

# In[2]:


paper_pmid = 23552365
paper_name = 'gonzalez_daran_2003' 


# In[3]:


datasets = pd.read_csv('extras/YeastPhenome_' + str(paper_pmid) + '_datasets_list.txt', sep='\t', header=None, names=['dataset_id', 'name'])


# In[4]:


datasets.set_index('dataset_id', inplace=True)


# # Load & process the data

# In[5]:


xl = pd.ExcelFile('raw_data/Butanol tolerance. Screening knockout collection.xlsx')
sheets = xl.sheet_names


# In[6]:


len(sheets)


# In[7]:


plate_list = []
for s in sheets:
    print('%s' % s)
    original_data = xl.parse(s)
    
    plate_num = int(s.split(' ')[0][2:])
    rows = np.arange(1,9)
    cols = np.arange(1,13)
        
    control_plate1 = pd.DataFrame(index=rows, columns=cols, data=original_data.iloc[5:13,2:14].values)
    test_plate1 = pd.DataFrame(index=rows, columns=cols, data=original_data.iloc[17:25,2:14].values)
    
    print(original_data.shape)
    
    if original_data.shape[1] > 14:
        control_plate2 = pd.DataFrame(index=rows, columns=cols, data=original_data.iloc[5:13,16:28].values)
        test_plate2 = pd.DataFrame(index=rows, columns=cols, data=original_data.iloc[17:25,16:28].values)
    else:
        control_plate2 = control_plate1.copy()
        test_plate2 = test_plate1.copy()
    
    control_plate1 = pd.melt(control_plate1.reset_index(), id_vars='index', var_name='col')
    control_plate1.columns = ['row','col','value']
    control_plate2 = pd.melt(control_plate2.reset_index(), id_vars='index', var_name='col')
    control_plate2.columns = ['row','col','value']
    
    control_plate = control_plate1.merge(control_plate2, on=['row','col'])
    control_plate['value'] = control_plate[['value_x','value_y']].mean(axis=1)
    
    test_plate1 = pd.melt(test_plate1.reset_index(), id_vars='index', var_name='col')
    test_plate1.columns = ['row','col','value']
    test_plate2 = pd.melt(test_plate2.reset_index(), id_vars='index', var_name='col')
    test_plate2.columns = ['row','col','value']
    
    test_plate = test_plate1.merge(test_plate2, on=['row','col'])
    test_plate['value'] = test_plate[['value_x','value_y']].mean(axis=1)
    
    plate = test_plate[['row','col','value']].merge(control_plate[['row','col','value']], on=['row','col'])
    plate['value'] = plate['value_x'] / plate['value_y']
    plate['num'] = plate_num
    
    plate_list.append(plate)


# In[8]:


original_data = pd.concat(plate_list, axis=0)


# In[9]:


original_data.shape


# In[10]:


# Load the platemap
plate_map = pd.read_excel('raw_data/Knockout collection map.xls', sheet_name='DATA')


# In[11]:


plate_map.shape


# In[12]:


plate_map = plate_map.loc[plate_map['Plate'].notnull() & plate_map['Row'].notnull() & plate_map['Col'].notnull(),:]


# In[13]:


plate_map['num'] = plate_map['Plate'].astype(int)


# In[14]:


plate_map['col'] = plate_map['Col'].astype(int)


# In[15]:


plate_map['row'] = plate_map['Row'].apply(lambda x: 'ABCDEFGH'.find(x)+1)


# In[16]:


original_data = original_data.merge(plate_map[['ORF name','num','col','row']], on=['num','row','col'])


# In[17]:


original_data.head()


# In[18]:


print('Original data dimensions: %d x %d' % (original_data.shape))


# In[19]:


original_data['orf'] = original_data['ORF name'].astype(str)


# In[20]:


# Eliminate all white spaces & capitalize
original_data['orf'] = clean_orf(original_data['orf'])


# In[21]:


original_data.loc[original_data['orf']=='YLR287-A','orf'] = 'YLR287C-A'


# In[22]:


# Translate to ORFs 
original_data['orf'] = translate_sc(original_data['orf'], to='orf')


# In[23]:


# Make sure everything translated ok
t = looks_like_orf(original_data['orf'])
print(original_data.loc[~t,])


# In[24]:


original_data['data'] = original_data['value']


# In[25]:


original_data.set_index('orf', inplace=True)


# In[26]:


original_data = original_data[['data']].copy()


# In[27]:


original_data = original_data.groupby(original_data.index).mean()


# In[28]:


original_data.shape


# # Prepare the final dataset

# In[29]:


data = original_data.copy()


# In[30]:


dataset_ids = [127]
datasets = datasets.reindex(index=dataset_ids)


# In[31]:


lst = [datasets.index.values, ['value']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data.columns = idx


# In[32]:


data.head()


# ## Subset to the genes currently in SGD

# In[33]:


genes = pd.read_csv(path_to_genes, sep='\t', index_col='id')
genes = genes.reset_index().set_index('systematic_name')
gene_ids = genes.reindex(index=data.index.values)['id'].values
num_missing = np.sum(np.isnan(gene_ids))
print('ORFs missing from SGD: %d' % num_missing)


# In[34]:


data['gene_id'] = gene_ids
data = data.loc[data['gene_id'].notnull()]
data['gene_id'] = data['gene_id'].astype(int)
data = data.reset_index().set_index(['gene_id','orf'])

data.head()


# # Normalize

# In[35]:


data_norm = normalize_phenotypic_scores(data, has_tested=True)


# In[36]:


# Assign proper column names
lst = [datasets.index.values, ['valuez']*datasets.shape[0]]
tuples = list(zip(*lst))
idx = pd.MultiIndex.from_tuples(tuples, names=['dataset_id','data_type'])
data_norm.columns = idx


# In[37]:


data_norm[data.isnull()] = np.nan
data_all = data.join(data_norm)

data_all.head()


# # Print out

# In[38]:


for f in ['value','valuez']:
    df = data_all.xs(f, level='data_type', axis=1).copy()
    df.columns = datasets['name'].values
    df = df.droplevel('gene_id', axis=0)
    df.to_csv(paper_name + '_' + f + '.txt', sep='\t')


# # Save to DB

# In[39]:


from IO.save_data_to_db3 import *


# In[40]:


save_data_to_db(data_all, paper_pmid)


# In[ ]:




