# -*- coding: utf-8 -*-
"""
Created on 13/07/2020

@author: yann-stanislas.barral@roche.com
"""

# Import the required packages
import sabs_pkpd
import matplotlib.pyplot as plt
import numpy as np
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
    return V / (96485 * 2.583592e-05) * 0.0001533576 - (Ki + Kss * 0.029411764705882353 + Nai + Nass * 0.029411764705882353 + 2*(tot_cai + tot_cass * 0.029411764705882353 + Cansr * 0.08117647059 + tot_cajsr * 0.007059) - Nao - 2 * Cao - Ko)


# Change Ki to 120 mM, Nai is set to its original value
init_cond = default[:]
init_cond[3] = 120 # Ki = 120 mM
init_cond[4] = 120 # Kss = 120 mM
init_cond[1] = 4 # Nai = 4 mM
init_cond[2] = 4 # Nass = 4 mM
s.set_default_state(init_cond)

# Set C0 value
ORd_maxi_G0 = G0_calc(Ki = init_cond[3], Kss = init_cond[4], Nai = init_cond[1], Nass = init_cond[2])
s.set_constant('membrane.c0', ORd_maxi_G0)

# Reset to initial conditions
s.reset()

# Run to limit cycle
s.pre(2000000)
a2 = s.run(1000, log_interval = 1)


# Change Ki to 150 mM, Nai is set to its original value
init_cond = default[:]
init_cond[3] = 152 # Ki = 150 mM
init_cond[4] = 152 # Kss = 150 mM
init_cond[1] = 16 # Nai = 4 mM
init_cond[2] = 16 # Nass = 4 mM
s.set_default_state(init_cond)

# Set C0 value
ORd_mini_G0 = G0_calc(Ki = init_cond[3], Kss = init_cond[4], Nai = init_cond[1], Nass = init_cond[2])
s.set_constant('membrane.c0', ORd_mini_G0)

# Reset to initial conditions
s.reset()

# Run to limit cycle
s.pre(2000000)
a3 = s.run(1000, log_interval = 1)


# In[For TT06]
# Define the function to compute C0 from the intracellular concentrations
def G0_calc(Ki, Nai, Cai = 0.00010578226860054304, Casr = 3.5556779165585235,
            Cass = 0.00021417287326980984, V=-84.9, Nao = 140, Ko = 5.4, Cao = 2):
    
    return V / 96.4853415 * 0.000185 / 0.016404 - (Ki + Nai + 2*(Cai + Cai*0.2 / (Cai + 0.001)) +
                                                 2*1.094/16.404*(Casr + Casr * 10 / (Casr + 0.3)) +
                                                 2*5.468e-2/16.404*(Cass + Cass * 0.4 / (Cass + 0.00025)) - Nao - Ko - 2 * Cao)

# Conservative change at start
init_cond = default1[:]
init_cond[5] = 120 # Ki = 120 mM
init_cond[4] = 4 # Nai = 4 mM
s1.set_default_state(init_cond)

# Set C0 value
TT06_maxi_G0 = G0_calc(Ki = init_cond[5], Nai = init_cond[4])
s1.set_constant('membrane.c0', TT06_maxi_G0)

# Reset to initial conditions
s1.reset()

# Run to limit cycle
s1.pre(2000000)
a = s1.run(1000, log_interval = 1)

# Conservative change at start
init_cond = default1[:]
init_cond[5] = 152 # Ki = 150 mM
init_cond[4] = 16 # Nai = 16 mM
s1.set_default_state(init_cond)

# Set C0 value
TT06_mini_G0 = G0_calc(Ki = init_cond[5], Nai = init_cond[4])
s1.set_constant('membrane.c0', TT06_mini_G0)

# Reset to initial conditions
s1.reset()

# Run to limit cycle
s1.pre(2000000)
a1 = s1.run(1000, log_interval = 1)


# In[Plot the resulting APs]
def place_caption_label(ax, label, loc='upper left', Fontsize=35):
    from matplotlib.offsetbox import AnchoredText
    at = AnchoredText(label, loc=loc, prop=dict(size=Fontsize), frameon=True, borderpad = 0)
    ax.add_artist(at)    
    return None

fig, ax = plt.subplots(1, 2, figsize = (7.5, 4))
size_labels = 15
size_legends = 12
color_low_Ki = 'orange'
color_high_Ki = 'darkgreen'
width = 3

# Plot the results for TTP06
index = 0
ax[index].plot(a['membrane.V'], label = '$\Gamma_0 = $' + str(np.round(TT06_maxi_G0, decimals = 1)) + ' mM', LineWidth = width, color = color_low_Ki)
ax[index].plot(a1['membrane.V'], label = '$\Gamma_0 = $' + str(np.round(TT06_mini_G0, decimals = 1)) + ' mM', LineWidth = width, color = color_high_Ki)
ax[index].set_xlabel('Time (ms)', Fontsize = size_labels)
ax[index].set_ylabel('Voltage (mV)', Fontsize = size_labels)
ax[index].legend(fontsize = size_legends)
place_caption_label(ax[index], 'A', loc = 'lower right', Fontsize = 25)

# Plot the results for ORd CiPA
index = 1
ax[index].plot(a2['membrane.V'], label = '$\Gamma_0 = $' + str(np.round(ORd_maxi_G0, decimals = 1)) + ' mM', LineWidth = width, color = color_low_Ki)
ax[index].plot(a3['membrane.V'], label = '$\Gamma_0 = $' + str(np.round(ORd_mini_G0, decimals = 1)) + ' mM', LineWidth = width, color = color_high_Ki)
ax[index].set_xlabel('Time (ms)', Fontsize = size_labels)
ax[index].set_ylabel('Voltage (mV)', Fontsize = size_labels)
ax[index].legend(fontsize = size_legends)
place_caption_label(ax[index], 'B', loc = 'lower right', Fontsize = 25)

plt.tight_layout()
plt.savefig('./Figures/APs_Gamma_0.png', dpi = 300)

