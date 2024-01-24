#!/usr/bin/env python
# coding: utf-8

# In[1]:


import dynamo as dyn
import os
import h5py
import anndata as ad
import seaborn as sns
import numpy as np
import pandas as pd
import scanpy as sc
import matplotlib.pyplot as plt
#dyn.get_all_dependencies_version()
os.chdir('/Users/bastiaanspanjaard/Documents/Projects/RNA_velocity/')


# # Own data

# In[2]:


# Load annotations
anno_1 = pd.read_csv('./Data/notch_pool_ann_transfer.csv', index_col = 'Barcode') # Use for b22 and b23
anno_2 = pd.read_csv('./Data/inj_pool_ann_transfer.csv', index_col = 'Barcode') # Use for b24 and b25
colors = pd.read_csv('./Data/mct_colors.csv', index_col = 'Celltype_major')


# P1407_NM1_S6 - control for Notch inhibition (DMSO treatment) - b22_tel  
# P1407_NM2_S7 - Notch inhibition (chemical perturbation) - b23_tel  
# P1407_NM3_S8 - control for injury - b24_tel  
# P1407_NM4_S9 - skull injury 3dpi - b25_tel  
# 
# All are obtained in the same way experimentally (6 hour labeling with 4SU prior to collecting the cells)

# In[4]:


# Load data
prefix = 'P1407_NM4_S9'
filename = './Data/' + prefix + '_estimates.h5ad'
adata = dyn.read_h5ad(filename)
adata.obs.index = 'B25.T_' + adata.obs.index + '-1'
adata.obs = adata.obs.join(anno_2)


# In[5]:


# Add cell types and select relevant ones
translations = pd.read_csv('./Data/rg_neu_translation_table.csv')
fourclass_dict = dict(zip(translations['Celltype_detailed'], translations['Grouped_4classes']))
sevenclass_dict = dict(zip(translations['Celltype_detailed'], translations['Grouped_7classes']))
adata.obs['Fourclass_types'] = adata.obs['Celltypes_detailed_pred'].replace(to_replace = fourclass_dict)
adata.obs['Sevenclass_types'] = adata.obs['Celltypes_detailed_pred'].replace(to_replace = sevenclass_dict)
adata_filter = adata[adata.obs['Sevenclass_types'].isin(translations['Grouped_7classes'].unique())]


# In[5]:


# Plot observed labeling rates
pe_plot = sns.histplot(adata_filter.obs, x = 'p_e').set_title(prefix)
pe_plot_fig = pe_plot.get_figure()
pe_plot_fig.savefig('./Images/' + prefix + 'p_e_histogram.png')


# In[6]:


# Plot estimated labeling rates
tc_plot = sns.histplot(adata_filter.obs, x = 'p_c_TC').set_title(prefix)
tc_plot_fig = tc_plot.get_figure()
tc_plot_fig.savefig('./Images/' + prefix + 'TC_histogram.png')


# In[7]:


# Create new adata object with correctly-named layers
#adata_li = adata_filter.copy()
adata_filter.X = adata_filter.layers['total']
adata_filter.layers['uu'] = adata_filter.layers['un_TC_est'].copy()
adata_filter.layers['ul'] = adata_filter.layers['ul_TC_est'].copy()
adata_filter.layers['sl'] = adata_filter.layers['sl_TC_est'].copy()
adata_filter.layers['su'] = adata_filter.layers['sn_TC_est'].copy()

for layer in list(adata_filter.layers.keys()):
    if layer not in ['uu', 'ul', 'sl', 'su']:
        del adata_filter.layers[layer]


# # Merge data

# In[ ]:


# See https://dynast-release.readthedocs.io/en/stable/technical_information.html#read-groups for explanation of the different layers.


# In[3]:


filename_1 = './Data/P1407_NM1_S6_estimates.h5ad'
adata_1 = dyn.read_h5ad(filename_1)
adata_1.obs.index = 'B22.T_' + adata_1.obs.index + '-1'
adata_1.obs = adata_1.obs.join(anno_1)
adata_1.X = adata_1.layers['total']
adata_1.layers['uu'] = adata_1.layers['un_TC_est'].copy()
adata_1.layers['ul'] = adata_1.layers['ul_TC_est'].copy()
adata_1.layers['sl'] = adata_1.layers['sl_TC_est'].copy()
adata_1.layers['su'] = adata_1.layers['sn_TC_est'].copy()

for layer in list(adata_1.layers.keys()):
    if layer not in ['uu', 'ul', 'sl', 'su']:
        del adata_1.layers[layer]


# In[4]:


filename_2 = './Data/P1407_NM2_S7_estimates.h5ad'
adata_2 = dyn.read_h5ad(filename_2)
adata_2.obs.index = 'B23.T_' + adata_2.obs.index + '-1'
adata_2.obs = adata_2.obs.join(anno_1)

adata_2.X = adata_2.layers['total']
adata_2.layers['uu'] = adata_2.layers['un_TC_est'].copy()
adata_2.layers['ul'] = adata_2.layers['ul_TC_est'].copy()
adata_2.layers['sl'] = adata_2.layers['sl_TC_est'].copy()
adata_2.layers['su'] = adata_2.layers['sn_TC_est'].copy()

for layer in list(adata_2.layers.keys()):
    if layer not in ['uu', 'ul', 'sl', 'su']:
        del adata_2.layers[layer]


# In[5]:


filename_3 = './Data/P1407_NM3_S8_estimates.h5ad'
adata_3 = dyn.read_h5ad(filename_3)
adata_3.obs.index = 'B24.T_' + adata_3.obs.index + '-1'
adata_3.obs = adata_3.obs.join(anno_2)

adata_3.X = adata_3.layers['total']
adata_3.layers['uu'] = adata_3.layers['un_TC_est'].copy()
adata_3.layers['ul'] = adata_3.layers['ul_TC_est'].copy()
adata_3.layers['sl'] = adata_3.layers['sl_TC_est'].copy()
adata_3.layers['su'] = adata_3.layers['sn_TC_est'].copy()

for layer in list(adata_3.layers.keys()):
    if layer not in ['uu', 'ul', 'sl', 'su']:
        del adata_3.layers[layer]


# In[6]:


filename_4 = './Data/P1407_NM4_S9_estimates.h5ad'
adata_4 = dyn.read_h5ad(filename_4)
adata_4.obs.index = 'B25.T_' + adata_4.obs.index + '-1'
adata_4.obs = adata_4.obs.join(anno_2)

adata_4.X = adata_4.layers['total']
adata_4.layers['uu'] = adata_4.layers['un_TC_est'].copy()
adata_4.layers['ul'] = adata_4.layers['ul_TC_est'].copy()
adata_4.layers['sl'] = adata_4.layers['sl_TC_est'].copy()
adata_4.layers['su'] = adata_4.layers['sn_TC_est'].copy()

for layer in list(adata_4.layers.keys()):
    if layer not in ['uu', 'ul', 'sl', 'su']:
        del adata_4.layers[layer]


# In[17]:


#adata = adata_1.concatenate(adata_3, adata_4) #adata_2, 
adata = adata_1.concatenate(adata_2)
prefix = 'Notch_inh_control'


# In[53]:


# Are the layers merged correctly? Looks like it.
merged_index = adata.var.index.intersection(adata_1.var.index)
merged_index_numbers = adata_1.var.index.get_indexer(merged_index.to_list())

ad_1_cells = [x + '-0' for x in adata_1.obs.index.to_list()]
ad_1_cells_index_numbers = adata.obs.index.get_indexer(ad_1_cells)
a = adata.layers['ul'][ad_1_cells_index_numbers, :]
b = adata_1.layers['ul'][:, merged_index_numbers]
(a!=b).nnz==0


# In[8]:


translations = pd.read_csv('./Data/rg_neu_translation_table.csv')
fourclass_dict = dict(zip(translations['Celltype_detailed'], translations['Grouped_4classes']))
sevenclass_dict = dict(zip(translations['Celltype_detailed'], translations['Grouped_7classes']))
adata.obs['Fourclass_types'] = adata.obs['Celltypes_detailed_pred'].replace(to_replace = fourclass_dict)
adata.obs['Sevenclass_types'] = adata.obs['Celltypes_detailed_pred'].replace(to_replace = sevenclass_dict)
adata_filter = adata[adata.obs['Sevenclass_types'].isin(translations['Grouped_7classes'].unique())]


# # Velocity based on labeling and intron reads

# In[9]:


adata_li = adata_filter.copy()
adata_li.obs['time'] = 6


# In[10]:


# Preprocess
dyn.pp.recipe_monocle(
    adata_li,
    tkey="time",
    experiment_type="one-shot")


# In[11]:


adata_li.uns["pp"]["has_splicing"] = True


# In[12]:


# Calculate dynamics
dyn.tl.dynamics(adata_li, group="time", one_shot_method="sci_fate", model="stochastic", re_smooth = True)


# In[13]:


sc.pp.neighbors(adata_li)


# In[14]:


sc.tl.umap(adata_li)


# In[15]:


dyn.tl.cell_velocities(adata_li, enforce = True)


# In[16]:


dyn.pl.streamline_plot(adata_li, vkey = 'velocity_T')


# In[18]:


adata_li_save = adata_li.copy()
adata_li_save.layers['velocity'] = adata_li_save.layers['velocity_T']
del adata_li_save.uns['velocyto_SVR']
del adata_li_save.uns['cell_phase_genes']
for layer in list(adata_li_save.layers.keys()):
    if layer not in ['new', 'total', 'spliced', 'unspliced', 'velocity']:
        del adata_li_save.layers[layer]
var_object_columns = [x for (x, b) in zip(list(adata_li_save.var), adata_li_save.var.dtypes == "object") if b]
adata_li_save.var[var_object_columns] = adata_li_save.var[var_object_columns].astype(float)
adata_li_save.write('./Data/' + prefix + '_intronSLAM_estimated.h5ad')
#adata_li_save.write('./Data/SLAM1_3_4_intronSLAM_estimated.h5ad')


# In[ ]:




