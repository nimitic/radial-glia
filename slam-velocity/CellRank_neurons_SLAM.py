#!/usr/bin/env python
# coding: utf-8

# In[37]:


from IPython.core.display import display, HTML
display(HTML("<style>.container { width:90% !important; }</style>"))
get_ipython().run_line_magic('matplotlib', 'inline')
import scanpy as sc
import numpy as np
import pandas as pd
import matplotlib.pyplot as pl
from matplotlib import rcParams
from sklearn import preprocessing
import seaborn as sns
import h5py
sc.set_figure_params(dpi_save = 300)
import scvelo as scv
from scipy import stats
scv.settings.set_figure_params('scvelo', dpi_save = 300)


# P1407_NM1_S6 - control for Notch inhibition (DMSO treatment) - b22_tel  
# P1407_NM2_S7 - Notch inhibition (chemical perturbation) - b23_tel  
# P1407_NM3_S8 - control for injury - b24_tel  
# P1407_NM4_S9 - skull injury 3dpi - b25_tel  
# 
# All are obtained in the same way experimentally (6 hour labeling with 4SU prior to collecting the cells)

# In[3]:


adata = sc.read_h5ad('../Data/Notch_inh_control_intronSLAM_estimated.h5ad') 
adata.var = adata.var.rename({'use_for_transition' : 'velocity_genes'}, axis = 1)
adata.layers['velocity'] = adata.layers['velocity'].toarray()
brain_area = "Telencephalon_Notch_inh_and_control" #"Telencephalon_Notch_inhibition_control"#"Telencephalon_4_SLAM"


# In[4]:


scv.pl.velocity_embedding_stream(adata, color = ['snap25a', 'fabp7a', 'tubb5', 'cd99l2',
                          'pcna', 'mki67', 'top2a', 'hmgn2', 'glula'], ncols = 3)#,
#          save = '_' + brain_area + '_intronSLAM_marker_genes.png')


# In[5]:


sc.tl.leiden(adata, resolution = 3)


# In[6]:


sc.set_figure_params(figsize = (8,6))
sc.pl.umap(adata, color = 'leiden', legend_loc='on data')
sc.set_figure_params(figsize = None)
scv.settings.set_figure_params('scvelo', dpi_save = 300)


# In[7]:


# DATASET Notch inhibition and inhibition control
# For resolution == 3
new_cluster_dict = {'0' : 'Neurons', '1' : 'Radial glia', '2' : 'Neurons',
                    '3' : 'Neurons newborn', '4' : 'Neurons', '5' : 'Neurons',
                    '6' : 'Radial glia', '7' : 'Neurons', '8' : 'Neurons',
                    '9' : 'Radial glia', '10' : 'Neurons', '11' : 'Neurons newborn',
                    '12' : 'Neurons', '13' : 'Radial glia', '14' : 'Neurons',
                    '15' : 'Neurons', '16' : 'Radial glia', '17' : 'Neurons',
                    '18' : 'Neurons', '19' : 'Neurons', '20': 'Neurons',
                    '21' : 'Neurons', '22' : 'Neurons', '23' : 'Neurons',
                    '24' : 'Neurons newborn', '25' : 'Proliferating cells', '26' : 'Neurons',
                    '27' : 'Radial glia', '28' : 'Neurons', '29' : 'Neurons newborn',
                    '30' : 'Radial glia', '31' : 'Radial glia', '32' : 'Proliferating cells',
                    '33' : 'Radial glia', '34' : 'Neurons', '35' : 'Radial glia',
                    '36' : 'Radial glia snap25a', '37' : 'Radial glia', '38' : 'Neurons',
                    '39' : 'Neurons', '40' : 'Neurons', '41' : 'Radial glia snap25a', 
                    '42' : 'Neurons', '43' : 'Neurons'
}

# DATASET P1407_NM2_S7/Telencephalon_Notch_inhibition
# For resolution == 3
# new_cluster_dict = {'0' : 'Neurons', '1' : 'Radial glia', '2' : 'Neurons',
#                     '3' : 'Neurons newborn', '4' : 'Neurons', '5' : 'Radial glia',
#                     '6' : 'Neurons', '7' : 'Neurons', '8' : 'Neurons',
#                     '9' : 'Neurons', '10' : 'Neurons', '11' : 'Neurons',
#                     '12' : 'Neurons', '13' : 'Neurons', '14' : 'Neurons',
#                     '15' : 'Neurons', '16' : 'Neurons', '17' : 'Neurons newborn',
#                     '18' : 'Radial glia', '19' : 'Radial glia', '20': 'Proliferating cells',
#                     '21' : 'Neurons', '22' : 'Radial glia', '23' : 'Neurons',
#                     '24' : 'Radial glia snap25a', '25' : 'Radial glia', '26' : 'Radial glia',
#                     '27' : 'Radial glia snap25a', '28' : 'Neurons newborn', '29' : 'Radial glia',
#                     '30' : 'Radial glia', '31' : 'Neurons', '32' : 'Neurons',
#                     '33' : 'Neurons', '34' : 'Neurons', '35' : 'Neurons',
#                     '36' : 'Neurons', '37' : 'Neurons', '38' : 'Neurons'
# }
# 0, 2, 4, 5, 7, 9, 10, 11, 12, 14, 15, 17, 18, 19, 22, 26, 28, 29, 32, 33, 34, 36 are Neurons
# 1, 6, 16, 20, 21, 25, 31 are Radial glia
# 3, 13(?), 23(?) are Neurons newborn
# 30 are proliferating cells
# 8 are radial glia her4++
# 24, 27, 35, 37 are radial glia snap25a

# DATASET P1407_NM1_S6/Telencephalon_Notch_control
# resolution == 3
# new_cluster_dict = {'0' : 'Neurons newborn', '1' : 'Neurons', '2' : 'Radial glia',
#                    '3' : 'Neurons', '4' : 'Neurons', '5' : 'Neurons',
#                    '6' : 'Neurons', '7' : 'Neurons', '8' : 'Neurons',
#                    '9' : 'Neurons', '10' : 'Neurons', '11' : 'Radial glia',
#                    '12' : 'Neurons', '13' : 'Radial glia', '14' : 'Radial glia',
#                    '15' : 'Neurons newborn', '16' : 'Radial glia', '17' : 'Neurons',
#                    '18' : 'Neurons newborn', '19' : 'Neurons', '20' : 'Neurons', 
#                     '21' : 'Proliferating cells', '22' : 'Neurons', '23' : 'Neurons',  # NB not sure whether 21 are really proliferating cells.
#                     '24' : 'Radial glia', '25' : 'Neurons', '26' : 'Radial glia', 
#                     '27' : 'Neurons', '28' : 'Neurons', '29' : 'Radial glia',
#                     '30' : 'Radial glia', '31' : 'Neurons', '32' : 'Radial glia',
#                     '33' : 'Radial glia', '34' : 'Radial glia', '35' : 'Neurons',
#                     '36' : 'Neurons', '37' : 'Radial glia snap25a', '38' : 'Proliferating cells', # NB not sure whether 38 are really proliferating cells.
#                     '39' : 'Neurons'}

adata.obs['Cell_type_mg'] = adata.obs['leiden'].replace(new_cluster_dict) # Mg = marker genes


# In[8]:


SLAM_markergene_colors = pd.DataFrame(data = ['#FAF9C8', '#F2EB3D', '#BABABA', '#F57F20','#E80C19'],
                                     index = ['Neurons', 'Neurons newborn', 'Radial glia', 'Proliferating cells', 'Radial glia snap25a'], columns = ['Color'])
SLAM_markergene_class_palette = SLAM_markergene_colors.reindex(adata.obs.Cell_type_mg.cat.categories)['Color'].to_list()


# In[9]:


scv.pl.velocity_embedding_stream(adata, color = 'Cell_type_mg', title = '', 
                                 legend_loc = 'right', #save = brain_area + '_intronSLAM_sc_umap_mg_types_res_3.png',
                                palette = SLAM_markergene_class_palette)


# In[10]:


# Plot separate datasets on combined umap embedding.
# P1407_NM1_S6 - control for Notch inhibition (DMSO treatment) - b22_tel  
# P1407_NM2_S7 - Notch inhibition (chemical perturbation) - b23_tel  
# P1407_NM3_S8 - control for injury - b24_tel  
# P1407_NM4_S9 - skull injury 3dpi - b25_tel  
# All are obtained in the same way experimentally (6 hour labeling with 4SU prior to collecting the cells)

sc.pl.umap(adata[adata.obs['Orig_ident'] == 'b23_tel'], color = 'Cell_type_mg', title = 'Inhibition', 
                                 legend_loc = 'right', save = 'Inhibition_on2_intronSLAM_sc_umap_mg_types_res_3.png')


# In[11]:


from cellrank.kernels import VelocityKernel
from cellrank.kernels import ConnectivityKernel
vk = VelocityKernel(adata).compute_transition_matrix()
ck = ConnectivityKernel(adata).compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
combined_kernel.compute_transition_matrix()


# In[12]:


from cellrank.estimators import GPCCA
g = GPCCA(combined_kernel) #.compute_schur(n_components=20)
cg = GPCCA(ck)


# In[13]:


cg.compute_schur(n_components=50)


# In[14]:


cg.plot_spectrum(n = 30, real_only = True)


# In[15]:


cg.compute_macrostates(n_states=17, cluster_key="Cell_type_mg")
cg.plot_macrostates(same_plot = False, ncols = 3, 
                    save = '_' + brain_area + '_macrostate_memberships.png')#,


# In[16]:


# Take macrostates from connectivity kernel and use them as terminal states for the connectivity + velocity kernel.
# This is the exact method used in the function set_terminal_states_from_macrostates.
memberships = cg.macrostates_memberships
states = cg.macrostates
colors = memberships[list(states.cat.categories)].colors
probs = (memberships.X / memberships.X.max(0)).max(1)
probs = pd.Series(probs, index=cg.adata.obs_names)
g._write_terminal_states(states, colors, probs, memberships, params = {})


# In[17]:


macrostate_memberships = pd.DataFrame(cg.macrostates_memberships.X, columns = cg.macrostates.cat.categories, index = cg.adata.obs.index)


# In[18]:


macrostate_memberships.to_csv('../Data/' + brain_area + '_macrostate_memberships.csv')


# In[19]:


g.compute_absorption_probabilities()
g.plot_absorption_probabilities(same_plot=False, ncols = 4,
                               save = '_' + brain_area + '_absorption_probabilities.png')


# In[20]:


g.absorption_probabilities


# In[21]:


g.absorption_probabilities.names


# In[22]:


absorption_probabilities = pd.DataFrame(g.absorption_probabilities.X, columns = cg.macrostates.cat.categories, index = g.adata.obs.index)
absorption_probabilities


# In[23]:


absorption_probabilities.to_csv('../Data/' + brain_area + '_absorption_probabilities.csv')


# In[24]:


# Merge macrostates and discriminate between no inhibition and Notch inhibition.
# Make dictionary of new states that contains the names of the states in them
state_dict = {'Radial glia' : ['Radial glia_1', 'Radial glia_2', 'Radial glia_3',
                              'Radial glia_4', 'Radial glia_5', 'Radial glia_6'],
             'Neurons 1' : ['Neurons_4', 'Neurons_5', 'Neurons_6', 'Neurons_7'],
              'Neurons 2' : ['Neurons_1'],
              'Neurons 3' : ['Neurons_2'],
              'Neurons 4' : ['Neurons_3'],
             'Neurons newborn' : ['Neurons newborn'],
             'Proliferating cells' : ['Proliferating cells_1', 'Proliferating cells_2'],
             'Radial glia snap25a+' : ['Radial glia snap25a']}

#['Radial glia_1', 'Radial glia_2', 'Radial glia snap25a',
#       'Proliferating cells_1', 'Radial glia_3', 'Radial glia_4',
#       'Radial glia_5', 'Neurons_1', 'Proliferating cells_2', 'Neurons_2',
#       'Neurons_3', 'Neurons_4', 'Neurons_5', 'Neurons newborn',
#       'Neurons_6', 'Radial glia_6', 'Neurons_7']


#state_dict = {'Radial glia' : ['Radial glia_1', 'Radial glia_2'],
#              'Neurons' : ['Neurons_1', 'Neurons_2', 'Neurons_3', 'Neurons_4', 'Neurons_5', 'Neurons_6', 'Neurons_7', 'Neurons_8'],
#             'Neurons newborn' : ['Neurons newborn_1', 'Neurons newborn_2'],
#             'Proliferating cells' : ['Proliferating cells'],
#             'Radial glia snap25a+' : ['Radial glia snap25a+_1', 'Radial glia snap25a+_2'],
#             'Radial glia id2b+' : ['Radial glia id2b+_1', 'Radial glia id2b+_2']}
#memberships = cg.macrostates_memberships
memberships = g.absorption_probabilities
sum_states_array = np.ndarray(shape = (memberships.base.shape[0], len(state_dict)))
sum_states_array.shape
ssa_loc = 0
# Loop over dictionary and create three things: an array with sums of previous states, a list of names and a list of colors (the first color in the base membership Lineage)
for new_state in state_dict:
#    print(new_state)
#    print(memberships[state_dict[new_state]].sum(axis = 1).tolist())
    sum_states_array[:,ssa_loc] = [x for [x] in memberships[state_dict[new_state]].sum(axis = 1).tolist()]
    ssa_loc+=1
sum_states_names = list(state_dict.keys())
color_indx = [memberships.names.tolist().index(x) for x in [state_dict[x][0] for x in state_dict.keys()]]
sum_states_colors = memberships.colors[color_indx]
#sum_states = Lineage(input_array = sum_states_array, names = sum_states_names, colors = sum_states_colors)
#sum_states.sum(axis = 1) # All sum to 1
new_absorption_probs = pd.DataFrame(data = sum_states_array, columns = sum_states_names, index = adata.obs.index)
new_absorption_probs = new_absorption_probs.merge(adata.obs['Sevenclass_types'], left_index = True, right_index = True)


# In[26]:


sum_states_colors


# In[27]:


dmso_index = adata.obs[adata.obs['Orig_ident'] == "b22_tel"].index
notch_index = adata.obs[adata.obs['Orig_ident'] == "b23_tel"].index


# In[28]:


new_absorption_probs_dmso_melt = new_absorption_probs.loc[dmso_index].melt(id_vars = 'Sevenclass_types', var_name = 'Terminal state', value_name = 'Transition')
new_absorption_probs_notch_melt = new_absorption_probs.loc[notch_index].melt(id_vars = 'Sevenclass_types', var_name = 'Terminal state', value_name = 'Transition')
new_absorption_probs_dmso_melt['Dummy'] = 'Hello'
new_absorption_probs_notch_melt['Dummy'] = 'Hello'


# In[29]:


g = sns.catplot(
    data=new_absorption_probs_dmso_melt, x = 'Dummy', y="Transition", col="Sevenclass_types", hue = 'Terminal state', col_wrap = 4,#x="Terminal_state",
    kind="bar", palette = sum_states_colors, height=3, aspect=1, legend=True,
)
g.set_axis_labels("", "Absorption probability")
g.set_xticklabels("")
g.set_titles("{col_name}")
#g.savefig("../Images/Absorption_probabilities_barplot_dmso.png")


# In[30]:


g = sns.catplot(
    data=new_absorption_probs_notch_melt, x = 'Dummy', y="Transition", col="Sevenclass_types", hue = 'Terminal state', col_wrap = 4,#x="Terminal_state",
    kind="bar", palette = sum_states_colors, height=3, aspect=1, legend=True,
)
g.set_axis_labels("", "Absorption probability")
g.set_xticklabels("")
g.set_titles("{col_name}")
#g.savefig("../Images/Absorption_probabilities_barplot_notch_inhib.png")


# In[46]:


nn_nn_notch = new_absorption_probs_notch_melt[(new_absorption_probs_notch_melt['Sevenclass_types'] == 'Neurons newborn') & (new_absorption_probs_notch_melt['Terminal state'] == 'Neurons newborn')]['Transition']
nn_nn_dmso = new_absorption_probs_dmso_melt[(new_absorption_probs_dmso_melt['Sevenclass_types'] == 'Neurons newborn') & (new_absorption_probs_dmso_melt['Terminal state'] == 'Neurons newborn')]['Transition']
stats.ttest_ind(nn_nn_notch, nn_nn_dmso, equal_var = False, alternative = 'greater')


# In[47]:


nn_nn_notch.mean()


# In[48]:


nn_nn_dmso.mean()


# In[49]:


pc_pc_notch = new_absorption_probs_notch_melt[(new_absorption_probs_notch_melt['Sevenclass_types'] == 'Proliferating cells') & (new_absorption_probs_notch_melt['Terminal state'] == 'Proliferating cells')]['Transition']
pc_pc_dmso = new_absorption_probs_dmso_melt[(new_absorption_probs_dmso_melt['Sevenclass_types'] == 'Proliferating cells') & (new_absorption_probs_dmso_melt['Terminal state'] == 'Proliferating cells')]['Transition']
stats.ttest_ind(pc_pc_notch, pc_pc_dmso, equal_var = False, alternative = 'greater')


# In[50]:


pc_pc_notch.mean()


# In[51]:


pc_pc_dmso.mean()


# In[ ]:




