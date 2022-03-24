# -*- coding: utf-8 -*-
"""
Created on 13/07/2020

@author: yann-stanislas.barral@roche.com
"""

# Import the required packages
import numpy as np
import sabs_pkpd
import matplotlib.pyplot as plt
import matplotlib.axes
from mpl_toolkits.mplot3d import Axes3D
import os

# Select the folder in which this repo is downloaded in the line below
os.chdir('The/location/of/the/root/folder/of/this/repo')


# In[Load the models]
s = sabs_pkpd.load_model.load_simulation_from_mmt('./Models/Ohara CiPA - analytical voltage.mmt')
s1 = sabs_pkpd.load_model.load_simulation_from_mmt('./Models/tentusscher_2006 - analytical voltage.mmt')

# Set the solver tolerance to fine tolerance
s.set_tolerance(abs_tol = 1e-08, rel_tol = 1e-06)
s1.set_tolerance(abs_tol = 1e-08, rel_tol = 1e-06)

# Retrieve the initial conditions from the models
default = s.state()
default1 = s1.state()


# In[For ORd CiPA]
# Define the function to compute C0 from the intracellular concentrations
def G0_calc(Ki = 144.65559, Kss = 144.65556, Nai = 7.2680045, Nass = 7.26809,
            Cai = 8.6e-5, Cansr = 1.61957, Cajsr = 1.571234014, Cass = 8.49e-5,
            V=-88, Nao = 140, Ko = 5.4, Cao = 1.8):
    tot_cai = Cai * (1 + 0.05 / (Cai + 0.00238) + 0.07/(Cai + 0.0005))
    tot_cass = Cass * (1 + 0.047 / (Cass + 0.00087) + 1.124/(Cass + 0.0087))
    tot_cajsr = Cajsr * (1 + 10 / (Cajsr + 0.8))
    G0 = V / (96485 * 2.583592e-05) * 0.0001533576 - (Ki + Kss * 0.029411764705882353 + Nai + Nass * 0.029411764705882353 + 2*(tot_cai + tot_cass * 0.029411764705882353 + Cansr * 0.08117647059 + tot_cajsr * 0.007059) - Nao - 2*Cao - Ko)

    print('G0 = ' + str(G0))
    return G0

# Change Ki to 120 mM, Nai is set to its original value
init_cond = default[:]
init_cond[3] = 120.5 # Ki = 120 mM
init_cond[4] = 120.5 # Kss = 120 mM
s.set_default_state(init_cond)

# Set C0 value
G0_low_ORd = G0_calc(Ki = init_cond[3], Kss = init_cond[4])
s.set_constant('membrane.c0', G0_low_ORd)

# Reset to initial conditions
s.reset()

# Run to limit cycle
s.pre(2000000)
a2 = s.run(1000, log_interval = 1)


# Change Ki to 150 mM, Nai is set to its original value
init_cond = default[:]
init_cond[3] = 149.6 # Ki = 150 mM
init_cond[4] = 149.6 # Kss = 150 mM
s.set_default_state(init_cond)

# Set C0 value
G0_high_ORd = G0_calc(Ki = init_cond[3], Kss = init_cond[4])
s.set_constant('membrane.c0', G0_high_ORd)

# Reset to initial conditions
s.reset()

# Run to limit cycle
s.pre(2000000)
a3 = s.run(1000, log_interval = 1)


# In[For TT06]
# Define the function to compute C0 from the intracellular concentrations
def G0_calc(Ki, Nai, Cai = 0.00010578226860054304, Casr = 3.5556779165585235,
            Cass = 0.00021417287326980984, V=-84.9, Nao = 140, Ko = 5.4, Cao = 2):
    
    G0 = V / 96.4853415 * 0.000185 / 0.016404 - (Ki + Nai + 2*(Cai + Cai*0.2 / (Cai + 0.001)) +
                                                 2*1.094/16.404*(Casr + Casr * 10 / (Casr + 0.3)) +
                                                 2*5.468e-2/16.404*(Cass + Cass * 0.4 / (Cass + 0.00025)) - Nao - Ko - 2 * Cao)

    print('G0 = ' + str(G0))
    return G0
    
# Conservative change at start
init_cond = default1[:]
init_cond[5] = 120 # Ki = 120 mM
s1.set_default_state(init_cond)

# Set C0 value
G0_low_TT = G0_calc(Ki = init_cond[5], Nai = init_cond[4])
s1.set_constant('membrane.c0', G0_low_TT)

# Reset to initial conditions
s1.reset()

# Run to limit cycle
s1.pre(2000000)
a = s1.run(1000, log_interval = 1)

# Conservative change at start
init_cond = default1[:]
init_cond[5] = 150 # Ki = 150 mM
s1.set_default_state(init_cond)

# Set C0 value
G0_high_TT = G0_calc(Ki = init_cond[5], Nai = init_cond[4])
s1.set_constant('membrane.c0', G0_high_TT)

# Reset to initial conditions
s1.reset()

# Run to limit cycle
s1.pre(2000000)
a1 = s1.run(1000, log_interval = 1)


# In[Read out the restitution portraits from the WebLab]

TT06_low_Ki = np.loadtxt('./WebLab results/Restitution TT06 low Ki.csv',
                         skiprows  = 0,
                         delimiter = ',')

TT06_high_Ki = np.loadtxt('./WebLab results/Restitution TT06 high Ki.csv',
                         skiprows  = 0,
                         delimiter = ',')

ORd_low_Ki = np.loadtxt('./WebLab results/Restitution ORd CiPA low Ki.csv',
                         skiprows  = 0,
                         delimiter = ',')

ORd_high_Ki = np.loadtxt('./WebLab results/Restitution ORd CiPA high Ki.csv',
                         skiprows  = 0,
                         delimiter = ',')


# In[Plot Results]
def place_caption_label(ax, label, loc='upper left', fontsize=35):
    from matplotlib.offsetbox import AnchoredText
    at = AnchoredText(label, loc=loc, prop=dict(size=fontsize), frameon=True, borderpad = 0)
    ax.add_artist(at)    
    return None

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
fig, ax = plt.subplots(2, 5, gridspec_kw={'width_ratios': [1, 1, 1, 1, 0.2]}, figsize = (18, 7))

size_labels = 15
size_legends = 15
size_ticks = 15
color_low_Ki = 'orange'
color_high_Ki = 'darkgreen'


# Plot the results for TTP06
label_low = '$\Gamma_0$ = ' + str(np.round(G0_low_TT, decimals = 1)) + ' mM'
label_high = '$\Gamma_0$ = ' + str(np.round(G0_high_TT, decimals = 1)) + ' mM'

ax[0, 0].plot(a1['ik1.IK1'], linewidth = 4, color = color_high_Ki, label = label_high)
ax[0, 0].plot(a['ik1.IK1'], linewidth = 4, color = color_low_Ki, label = label_low)
ax[0, 0].set_xlabel('Time (ms)', fontsize = size_labels)
ax[0, 0].set_ylabel('$I_{K1}$ (pA/pF)', fontsize = size_labels)
ax[0, 0].set_xticks([0, 250, 500, 750, 1000])
ax[0, 0].set_xticklabels([0, 250, 500, 750, 1000])
ax[0, 0].tick_params(labelsize = size_ticks)
place_caption_label(ax[0, 0], 'A', loc='upper right', fontsize=25)

ax[0, 1].plot(a1['inak.INaK'], linewidth = 4, color = color_high_Ki)
ax[0, 1].plot(a['inak.INaK'], linewidth = 4, color = color_low_Ki)
ax[0, 1].set_xlabel('Time (ms)', fontsize = size_labels)
ax[0, 1].set_ylabel('$I_{NaK}$ (pA/pF)', fontsize = size_labels)
ax[0, 1].set_xticks([0, 250, 500, 750, 1000])
ax[0, 1].set_xticklabels([0, 250, 500, 750, 1000])
ax[0, 1].tick_params(labelsize = size_ticks)
place_caption_label(ax[0, 1], 'B', loc='upper right', fontsize=25)

ax[0, 2].plot(a1['membrane.V'], linewidth = 4, color = color_high_Ki, label = label_high)
ax[0, 2].plot(a['membrane.V'], linewidth = 4, color = color_low_Ki, label = label_low)
ax[0, 2].set_xlabel('Time (ms)', fontsize = size_labels)
ax[0, 2].set_ylabel('Voltage (mV)', fontsize = size_labels)
ax[0, 2].set_xticks([0, 250, 500, 750, 1000])
ax[0, 2].set_xticklabels([0, 250, 500, 750, 1000])
ax[0, 2].tick_params(labelsize = size_ticks)
place_caption_label(ax[0, 2], 'C', loc='upper right', fontsize=25)

# Plot restitution portrait for TTP06
ax[0, 3].plot(TT06_low_Ki[:, 0], TT06_low_Ki[:, 1], color = 'orange', 
         linewidth = 4)
ax[0, 3].plot(TT06_low_Ki[:, 0], TT06_low_Ki[:, 2], color = 'orange', linewidth = 4,
         label = label_low)

ax[0, 3].plot(TT06_high_Ki[:, 0], TT06_high_Ki[:, 1], color = 'darkgreen', linewidth = 4)
ax[0, 3].plot(TT06_high_Ki[:, 0], TT06_high_Ki[:, 2], color = 'darkgreen',
         linewidth = 4, label = label_high)

# Set Figure labels and legend
ax[0, 3].set_ylabel('$APD_{90} (ms)$', fontsize = size_labels)
ax[0, 3].set_xlabel('Pacing period (ms)', fontsize = size_labels)
ax[0, 3].set_xticks([250, 500, 1000, 1500, 2000])
ax[0, 3].set_xticklabels([250, '  500', 1000, 1500, 2000])
ax[0, 3].tick_params(labelsize = size_ticks)
place_caption_label(ax[0, 3], 'D', loc='upper right', fontsize=25)

# Add an inset to zoom into the short pacing periods
axins1 = inset_axes(ax[0, 3], width="50%", height="50%", loc=4, borderpad=2.5)

axins1.plot(TT06_high_Ki[:, 0], TT06_high_Ki[:, 1], color = 'darkgreen', linewidth = 4)
axins1.plot(TT06_high_Ki[:, 0], TT06_high_Ki[:, 2], color = 'darkgreen',
         linewidth = 4, label = label_high)
axins1.plot(TT06_low_Ki[:, 0], TT06_low_Ki[:, 1], color = 'orange', 
         linewidth = 4)
axins1.plot(TT06_low_Ki[:, 0], TT06_low_Ki[:, 2], color = 'orange', linewidth = 4,
         label = label_low)

# set up the inset
axins1.set_xticks([250, 500, 1000, 1500, 2000])
axins1.tick_params(labelsize = 17)
axins1.set_xlim(250, 500)
axins1.set_ylim(200, 275)

# Add the legend
ax[0, 4].axis('off')
ax[0, 4].plot([], [], color = 'darkgreen', linewidth = 4, label = label_high)
ax[0, 4].plot([], [], color = 'orange', linewidth = 4, label = label_low)
ax[0, 4].legend(fontsize = 17, loc = 'center left', bbox_to_anchor=(-2,0.5))



# Plot the results for ORd CiPA
label_low = '$\Gamma_0$ = ' + str(np.round(G0_low_ORd, decimals = 1)) + ' mM'
label_high = '$\Gamma_0$ = ' + str(np.round(G0_high_ORd, decimals = 1)) + ' mM'

ax[1, 0].plot(a3['IK1.IK1'], linewidth = 4, color = color_high_Ki)
ax[1, 0].plot(a2['IK1.IK1'], linewidth = 4, color = color_low_Ki)
ax[1, 0].set_xlabel('Time (ms)', fontsize = size_labels)
ax[1, 0].set_ylabel('$I_{K1}$ (pA/pF)', fontsize = size_labels)
ax[1, 0].set_xticks([0, 250, 500, 750, 1000])
ax[1, 0].set_xticklabels([0, 250, 500, 750, 1000])
ax[1, 0].tick_params(labelsize = size_ticks)
place_caption_label(ax[1, 0], 'E', loc='upper right', fontsize=25)

ax[1, 1].plot(a3['INaK.INaK'], linewidth = 4, color = color_high_Ki)
ax[1, 1].plot(a2['INaK.INaK'], linewidth = 4, color = color_low_Ki)
ax[1, 1].set_xlabel('Time (ms)', fontsize = size_labels)
ax[1, 1].set_ylabel('$I_{NaK}$ (pA/pF)', fontsize = size_labels)
ax[1, 1].set_xticks([0, 250, 500, 750, 1000])
ax[1, 1].set_xticklabels([0, 250, 500, 750, 1000])
ax[1, 1].tick_params(labelsize = size_ticks)
place_caption_label(ax[1, 1], 'F', loc='upper right', fontsize=25)

ax[1, 2].plot(a3['membrane.V'], linewidth = 4, color = color_high_Ki)
ax[1, 2].plot(a2['membrane.V'], linewidth = 4, color = color_low_Ki)
ax[1, 2].set_xlabel('Time (ms)', fontsize = size_labels)
ax[1, 2].set_ylabel('Voltage (mV)', fontsize = size_labels)
ax[1, 2].set_xticks([0, 250, 500, 750, 1000])
ax[1, 2].set_xticklabels([0, 250, 500, 750, 1000])
ax[1, 2].tick_params(labelsize = size_ticks)
place_caption_label(ax[1, 2], 'G', loc='upper right', fontsize=25)


# Plot restitution portrait for ORd CiPA
ax[1, 3].plot(ORd_low_Ki[:, 0], ORd_low_Ki[:, 1], color = 'orange', 
         linewidth = 4)
ax[1, 3].plot(ORd_low_Ki[:, 0], ORd_low_Ki[:, 2], color = 'orange', linewidth = 4,
         label = label_low)

ax[1, 3].plot(ORd_high_Ki[:, 0], ORd_high_Ki[:, 1], color = 'darkgreen', linewidth = 4)
ax[1, 3].plot(ORd_high_Ki[:, 0], ORd_high_Ki[:, 2], color = 'darkgreen',
         linewidth = 4, label = label_high)

# Set up the labels and legend
ax[1, 3].set_ylabel('$APD_{90} (ms)$', fontsize = size_labels)
ax[1, 3].set_xlabel('Pacing period (ms)', fontsize = size_labels)
ax[1, 3].set_xticks([250, 500, 1000, 1500, 2000])
ax[1, 3].set_xticklabels([250, '  500', 1000, 1500, 2000])
ax[1, 3].tick_params(labelsize = size_ticks)
#ax[1, 3].legend(fontsize = 25, loc = 'upper left')
place_caption_label(ax[1, 3], 'H', loc='upper right', fontsize=25)

# Add an inset to zoom into the short pacing periods
axins2 = inset_axes(ax[1, 3], width="50%", height="50%", loc=4, borderpad=2.5)

axins2.plot(ORd_low_Ki[:, 0], ORd_low_Ki[:, 1], color = 'orange', 
         linewidth = 4)
axins2.plot(ORd_low_Ki[:, 0], ORd_low_Ki[:, 2], color = 'orange', linewidth = 4,
         label = label_low)

axins2.plot(ORd_high_Ki[:, 0], ORd_high_Ki[:, 1], color = 'darkgreen', linewidth = 4)
axins2.plot(ORd_high_Ki[:, 0], ORd_high_Ki[:, 2], color = 'darkgreen',
         linewidth = 4, label = label_high)

# set up the inset
axins2.set_xticks([250, 500, 1000, 1500, 2000])
axins2.tick_params(labelsize = 17)
axins2.set_xlim(250, 500)
axins2.set_ylim(180, 220)

# Add the legend
ax[1, 4].axis('off')
ax[1, 4].plot([], [], color = 'darkgreen', linewidth = 4, label = label_high)
ax[1, 4].plot([], [], color = 'orange', linewidth = 4, label = label_low)
ax[1, 4].legend(fontsize = 17, loc = 'center left', bbox_to_anchor=(-2,0.5))

plt.tight_layout()
plt.savefig('./Figures/outputs_and_Gamma_0.png', dpi = 300)

