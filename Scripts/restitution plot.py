# -*- coding: utf-8 -*-
"""
Created on 13/07/2020

@author: yann-stanislas.barral@roche.com
"""

# Import the required packages
import numpy as np
import matplotlib.pyplot as plt

# Load the results output from the WebLab
folder = 'C:/Users/barraly/Documents/PhD/Sensitivity to initial conditions/initial_conditions_study/'
TT06_low_Ki = np.loadtxt(folder + 'WebLab results/Restitution TT06 low Ki.csv',
                         skiprows  = 0,
                         delimiter = ',')

TT06_high_Ki = np.loadtxt(folder + 'WebLab results/Restitution TT06 high Ki.csv',
                         skiprows  = 0,
                         delimiter = ',')

ORd_low_Ki = np.loadtxt(folder + 'WebLab results/Restitution ORd CiPA low Ki.csv',
                         skiprows  = 0,
                         delimiter = ',')

ORd_high_Ki = np.loadtxt(folder + 'WebLab results/Restitution ORd CiPA high Ki.csv',
                         skiprows  = 0,
                         delimiter = ',')


# In[Plot Results]
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
fig, ax = plt.subplots(2, 1, figsize=[7, 14])

# Plot for TTP06
ax[0].plot(TT06_low_Ki[:, 0], TT06_low_Ki[:, 1], color = 'orange', 
         LineWidth = 4)
ax[0].plot(TT06_low_Ki[:, 0], TT06_low_Ki[:, 2], color = 'orange', LineWidth = 4,
         label = '$\Gamma_0=17.6$ mM')

ax[0].plot(TT06_high_Ki[:, 0], TT06_high_Ki[:, 1], color = 'darkgreen', LineWidth = 4)
ax[0].plot(TT06_high_Ki[:, 0], TT06_high_Ki[:, 2], color = 'darkgreen',
         LineWidth = 4, label = '$\Gamma_0=-13.3$ mM')

# Set Figure labels and legend
ax[0].set_ylabel('$APD_{90} (ms)$', Fontsize = 30)
ax[0].set_xlabel('Pacing period (ms)', Fontsize = 30)
ax[0].set_xticks([250, 500, 1000, 1500, 2000])
ax[0].set_xticklabels([250, '  500', 1000, 1500, 2000])
ax[0].tick_params(labelsize = 28)
ax[0].legend(fontsize = 25, loc = 'upper left')

# Add an inset to zoom into the short pacing periods
axins1 = inset_axes(ax[0], width="50%", height="50%", loc=4, borderpad=2.5)

axins1.plot(TT06_high_Ki[:, 0], TT06_high_Ki[:, 1], color = 'darkgreen', LineWidth = 4)
axins1.plot(TT06_high_Ki[:, 0], TT06_high_Ki[:, 2], color = 'darkgreen',
         LineWidth = 4, label = '$\Gamma_0=-13.3$ mM')
axins1.plot(TT06_low_Ki[:, 0], TT06_low_Ki[:, 1], color = 'orange', 
         LineWidth = 4)
axins1.plot(TT06_low_Ki[:, 0], TT06_low_Ki[:, 2], color = 'orange', LineWidth = 4,
         label = '$\Gamma_0=17.6$ mM')

# set up the inset
axins1.set_xticks([250, 500, 1000, 1500, 2000])
axins1.tick_params(labelsize = 17)
axins1.set_xlim(250, 500)
axins1.set_ylim(200, 275)


# Plot for ORd CiPA
ax[1].plot(ORd_low_Ki[:, 0], ORd_low_Ki[:, 1], color = 'orange', 
         LineWidth = 4)
ax[1].plot(ORd_low_Ki[:, 0], ORd_low_Ki[:, 2], color = 'orange', LineWidth = 4,
         label = '$\Gamma_0=17.6$ mM')

ax[1].plot(ORd_high_Ki[:, 0], ORd_high_Ki[:, 1], color = 'darkgreen', LineWidth = 4)
ax[1].plot(ORd_high_Ki[:, 0], ORd_high_Ki[:, 2], color = 'darkgreen',
         LineWidth = 4, label = '$\Gamma_0=-13.3$ mM')

# Set up the labels and legend
ax[1].set_ylabel('$APD_{90} (ms)$', Fontsize = 28)
ax[1].set_xlabel('Pacing period (ms)', Fontsize = 28)
ax[1].set_xticks([250, 500, 1000, 1500, 2000])
ax[1].set_xticklabels([250, '  500', 1000, 1500, 2000])
ax[1].tick_params(labelsize = 28)
ax[1].legend(fontsize = 25, loc = 'upper left')

# Add an inset to zoom into the short pacing periods
axins2 = inset_axes(ax[1], width="50%", height="50%", loc=4, borderpad=2.5)

axins2.plot(ORd_low_Ki[:, 0], ORd_low_Ki[:, 1], color = 'orange', 
         LineWidth = 4)
axins2.plot(ORd_low_Ki[:, 0], ORd_low_Ki[:, 2], color = 'orange', LineWidth = 4,
         label = '$\Gamma_0=17.6$ mM')

axins2.plot(ORd_high_Ki[:, 0], ORd_high_Ki[:, 1], color = 'darkgreen', LineWidth = 4)
axins2.plot(ORd_high_Ki[:, 0], ORd_high_Ki[:, 2], color = 'darkgreen',
         LineWidth = 4, label = '$\Gamma_0=-13.3$ mM')

# set up the inset
axins2.set_xticks([250, 500, 1000, 1500, 2000])
axins2.tick_params(labelsize = 17)
axins2.set_xlim(250, 500)
axins2.set_ylim(180, 220)

plt.tight_layout()

