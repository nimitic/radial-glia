#!/usr/bin/env python
# coding: utf-8

# In[1]:


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
import cellrank as cr
scv.settings.set_figure_params('scvelo', dpi_save = 300)


# In[2]:


translations = pd.read_csv('../Data/rg_neu_translation_table.csv')
fourclass_dict = dict(zip(translations['Celltype_detailed'], translations['Grouped_4classes']))
sevenclass_dict = dict(zip(translations['Celltype_detailed'], translations['Grouped_7classes']))


# In[3]:


adata = sc.read_h5ad('../Data/tel_rg_neu_scv')
brain_area = "Telencephalon"


# In[4]:


adata.obs['Fourclass_types'] = adata.obs['celltype_detailed'].replace(to_replace = fourclass_dict)
adata.obs['Sevenclass_types'] = adata.obs['celltype_detailed'].replace(to_replace = sevenclass_dict)


# In[5]:


sevenclass_palette = ['#FAF9C8', '#F2EB3D', '#F57F20', '#984F9F', '#327EBA', '#BABABA', '#E80C19']


# In[6]:


sc.pl.umap(adata, color = 'Sevenclass_types', save = '_' + brain_area + '_sevenclass_types.png',
          palette = sevenclass_palette)


# In[7]:


scv.pl.velocity_embedding_stream(adata, color = ['snap25a', 'fabp7a', 'tubb5', 'cd99l2',
                          'pcna', 'mki67', 'top2a', 'hmgn2', 'her4.1', 'her4.4', 'glula'], ncols = 3,
          save = '_' + brain_area + '_marker_genes.png')


# In[8]:


scv.pl.velocity_embedding_stream(adata, color = 'Sevenclass_types', title = '', 
                                 legend_loc = 'right',
                                save = brain_area + '_sevenclass_types.png')


# In[7]:


from cellrank.kernels import VelocityKernel
from cellrank.kernels import ConnectivityKernel
vk = VelocityKernel(adata).compute_transition_matrix()
ck = ConnectivityKernel(adata).compute_transition_matrix()
combined_kernel = 0.8 * vk + 0.2 * ck
combined_kernel.compute_transition_matrix()


# In[8]:


from cellrank.estimators import GPCCA
g = GPCCA(combined_kernel) #.compute_schur(n_components=20)
cg = GPCCA(ck)


# In[9]:


cg.compute_schur(n_components=50)


# In[10]:


#g.plot_spectrum(n = 30, real_only = True)
cg.plot_spectrum(n = 30, real_only = True)


# In[11]:


# g.compute_macrostates(n_states=8, cluster_key="Sevenclass_types")
# g.plot_macrostates(same_plot = False, ncols = 3)
cg.compute_macrostates(n_states=17, cluster_key="Sevenclass_types")
cg.plot_macrostates(same_plot = False, ncols = 3)#,
                   #save = brain_area + '_seventypes_conn_macrostates.png')


# In[12]:


# Use collapsed macrostates from connectivity kernel and use them as terminal states for the connectivity + velocity kernel.
# This is the exact method used in the function set_terminal_states_from_macrostates.
memberships = cg.macrostates_memberships
states = cg.macrostates
colors = memberships[list(states.cat.categories)].colors
probs = (memberships.X / memberships.X.max(0)).max(1)
probs = pd.Series(probs, index=cg.adata.obs_names)
g._write_terminal_states(states, colors, probs, memberships, params = {})


# In[13]:


macrostate_memberships = pd.DataFrame(cg.macrostates_memberships.X, columns = cg.macrostates.cat.categories, index = cg.adata.obs.index)
macrostate_memberships.to_csv('../Data/' + brain_area + '_macrostate_memberships.csv')


# In[14]:


g.compute_absorption_probabilities()
g.plot_absorption_probabilities(same_plot=False, ncols = 4)#,
                               #save = '_' + brain_area + '_absorption_probabilities.png')


# In[16]:


absorption_probabilities = pd.DataFrame(g.absorption_probabilities.X, columns = cg.macrostates.cat.categories, index = g.adata.obs.index)
#absorption_probabilities.to_csv('../Data/' + brain_area + '_absorption_probabilities.csv')


# In[17]:


# Make dictionary of new states that contains the names of the states in them
state_dict = {'Radial glia' : ['Radial glia_1', 'Radial glia_2'],
              'Neurons' : ['Neurons_1', 'Neurons_2', 'Neurons_3', 'Neurons_4', 'Neurons_5', 'Neurons_6', 'Neurons_7', 'Neurons_8'],
             'Neurons newborn' : ['Neurons newborn_1', 'Neurons newborn_2'],
             'Proliferating cells' : ['Proliferating cells'],
             'Radial glia snap25a+' : ['Radial glia snap25a+_1', 'Radial glia snap25a+_2'],
             'Radial glia id2b+' : ['Radial glia id2b+_1', 'Radial glia id2b+_2']}
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


# In[19]:


new_absorption_probs_melt = new_absorption_probs.melt(id_vars = 'Sevenclass_types', var_name = 'Terminal state', value_name = 'Transition')


# In[20]:


new_absorption_probs_melt['Dummy'] = 'Hello'


# In[21]:


absfig = sns.catplot(
    data=new_absorption_probs_melt, x = 'Dummy', y="Transition", col="Sevenclass_types", hue = 'Terminal state', col_wrap = 4,#x="Terminal_state",
    kind="bar", palette = sum_states_colors, height=3, aspect=1, legend=True,
)
absfig.set_axis_labels("", "Absorption probability")
absfig.set_xticklabels("")
absfig.set_titles("{col_name}")
#absfig.savefig("../Images/Absorption_probabilities_barplot_six_terminal_states.png")


# In[ ]:




