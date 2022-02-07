# -*- coding: utf-8 -*-
"""
Created on Thu Sep  3 14:44:35 2020

@author: barraly
"""

import sabs_pkpd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.colors as clrs
import time
import os

# Select the folder in which this repo is downloaded in the line below
os.chdir('The/location/of/the/root/folder/of/this/repo')


# In[Load the models for TTP06]
# TTP06 with analytical voltage expression
analytic_TT06 = sabs_pkpd.load_model.load_simulation_from_mmt('./Models/tentusscher_2006 - analytical voltage.mmt')

# TTP06 with derivative voltage expression
derivative_TT06 = sabs_pkpd.load_model.load_simulation_from_mmt('./Models/tentusscher_2006.mmt')

readout_TT06 = 'potassium.K_i'


# In[Load the models for ORd]
# ORd with analytical voltage expression
analytic_ORd = sabs_pkpd.load_model.load_simulation_from_mmt('./Models/Ohara CiPA - analytical voltage.mmt')

# ORd with derivative voltage expression
derivative_ORd = sabs_pkpd.load_model.load_simulation_from_mmt('./Models/OHara CiPA.mmt')

readout_ORd = 'intracellular_ions.ki'


# In[Try different solver tolerances and plot evolution of state]
solver_tolerances = [1e-09, 1e-08, 1e-07, 1e-06, 1e-05]

# Run the model to paced limit cycle with the various solver tolerances
for tol in solver_tolerances:
    # Run the analytical model to limit cycle
    analytic.reset()
    analytic.set_tolerance(abs_tol = 1e-08, rel_tol = tol)
    an = analytic.run(3000000, log_interval = 1000)

    # Run the derivative model to limit cycle
    derivative.reset()
    derivative.set_tolerance(abs_tol = 1e-07, rel_tol = tol)
    der = derivative.run(3000000, log_interval = 1000)
    
    # Plot the convergence to limit cycle and stability after the 2000th pace
    plt.figure(figsize = (15, 10))
    plt.subplot(1, 2, 1)
    plt.plot(np.linspace(0, 3000, 3000), an[readout], label = 'Analytical')
    plt.plot(np.linspace(0, 3000, 3000), der[readout], label = 'Derivative')
    plt.legend(fontsize = 15)
    plt.xlabel('# Paces', fontsize = 17)
    plt.ylabel('$K_i$ (mM)', fontsize = 17)
    plt.title('Solver tolerance : ' + str(tol), fontsize = 24)
    
    plt.subplot(1, 2, 2)
    plt.plot(np.linspace(2000, 3000, 1000), an[readout][2000:], label = 'Analytical')
    plt.plot(np.linspace(2000, 3000, 1000), der[readout][2000:], label = 'Derivative')
    plt.legend(fontsize = 15)
    plt.xlabel('# Paces', fontsize = 17)
    plt.ylabel('$K_i$ (mM)', fontsize = 17)
    plt.title('Solver tolerance : ' + str(tol), fontsize = 24)
    
    plt.tight_layout()
    

# In[Try different solver tolerances]
def maps(analytic, derivative, readout):
    solver_tolerances = [1e-09, 1e-08, 1e-07, 1e-06, 1e-05]
    stability_analytic = np.zeros((len(solver_tolerances), len(solver_tolerances)))
    stability_derivative = np.zeros((len(solver_tolerances), len(solver_tolerances)))
    durations_of_sim_an = np.zeros((len(solver_tolerances), len(solver_tolerances)))
    durations_of_sim_der = np.zeros((len(solver_tolerances), len(solver_tolerances)))
    
    # Run the model to paced limit cycle with the various solver tolerances
    for i, relative in enumerate(solver_tolerances):
        for j, absolute in enumerate(solver_tolerances):
            # Run the analytical model to limit cycle
            analytic.reset()
            analytic.set_tolerance(abs_tol = absolute, rel_tol = relative)
            beginning = time.time()
            an = analytic.run(3000000, log_interval = 1000)
            end = time.time()
            potassium_an = np.array(an[readout][2000:])
            
            print('It took ' + 
                  str(end - beginning) +
                  ' s to run the analytical model with ' +
                  str(absolute) + ' absolute and ' +
                  str(relative) + ' relative solver tolerance.\n')
            
            durations_of_sim_an[i, j] = end-beginning
            
            # Compute the stabiliy defined as mean distance to the value after 2k paces
            stability_analytic[i, j] = np.sum(abs(potassium_an - an[readout][2000])) / 1000
            
            # Run the derivative model to limit cycle
            derivative.reset()
            derivative.set_tolerance(abs_tol = absolute, rel_tol = relative)
            beginning = time.time()
            der = derivative.run(3000000, log_interval = 1000)
            end = time.time()
            
            potassium_der = np.array(der[readout][2000:])
            
            print('It took ' + 
                  str(end - beginning) +
                  ' s to run the derivative model with ' +
                  str(absolute) + ' absolute and ' +
                  str(relative) + ' relative solver tolerance.\n')
            
            durations_of_sim_der[i, j] = end-beginning
            
            # Compute the stabiliy defined as mean distance to the value after 2k paces
            stability_derivative[i, j] = np.sum(abs(potassium_der - der[readout][2000])) / 1000
    
    return stability_analytic, stability_derivative, durations_of_sim_an, durations_of_sim_der

# Compute the maps
stability_analytic_TT06, stability_derivative_TT06, durations_of_sim_an_TT06, durations_of_sim_der_TT06 = maps(analytic_TT06, derivative_TT06, readout_TT06)
stability_analytic_ORd, stability_derivative_ORd, durations_of_sim_an_ORd, durations_of_sim_der_ORd= maps(analytic_ORd, derivative_ORd, readout_ORd)

               
# In[ Plot the convergence to limit cycle and stability after the 2000th pace]
from matplotlib.gridspec import GridSpec
def place_caption_label(ax, label, loc='upper left', fontsize=35):
    from matplotlib.offsetbox import AnchoredText
    at = AnchoredText(label, loc=loc, prop=dict(size=fontsize), frameon=True, borderpad = 0)
    ax.add_artist(at)
    return None

fig = plt.figure(constrained_layout=False, figsize = (24, 18))

gs = GridSpec(3, 4, figure=fig, width_ratios = [0.4, 1, 1, 0.1], height_ratios = [0.2, 1, 1],
              hspace = 0.1, wspace = 0.1)

# For the maps
ax1 = fig.add_subplot(gs[1, 1])
plt.imshow(np.log10(stability_analytic_TT06))
ax1.set_ylim([-0.5, 4.5])
plt.clim(-8.5, -2.5)


ax2 = fig.add_subplot(gs[1, 2])
plt.imshow(np.log10(stability_derivative_TT06))
ax2.set_ylim([-0.5, 4.5])
plt.clim(-8.5, -2.5)

ax3 = fig.add_subplot(gs[2, 1])
plt.imshow(np.log10(stability_analytic_ORd))
ax3.set_ylim([-0.5, 4.5])
plt.clim(-8.5, -2.5)

ax4 = fig.add_subplot(gs[2, 2])
image = plt.imshow(np.log10(stability_derivative_ORd))
ax4.set_ylim([-0.5, 4.5])
plt.clim(-8.5, -2.5)

# Label axes
ticksize = 22
labelsize = 28
ax1.set_ylabel('Relative solver tolerance', fontsize = labelsize)
ax3.set_ylabel('Relative solver tolerance', fontsize = labelsize)
ax1.set_yticks([0, 1, 2, 3, 4])
ax1.set_yticklabels(['$10^{-9}$', '$10^{-8}$', '$10^{-7}$', '$10^{-6}$', '$10^{-5}$'])
ax3.set_yticks([0, 1, 2, 3, 4])
ax3.set_yticklabels(['$10^{-9}$', '$10^{-8}$', '$10^{-7}$', '$10^{-6}$', '$10^{-5}$'])
ax2.set_yticks([])
ax4.set_yticks([])

ax3.set_xlabel('Absolute solver tolerance', fontsize = labelsize)
ax4.set_xlabel('Absolute solver tolerance', fontsize = labelsize)
ax3.set_xticks([0, 1, 2, 3, 4])
ax3.set_xticklabels(['$10^{-9}$', '$10^{-8}$', '$10^{-7}$', '$10^{-6}$', '$10^{-5}$'])
ax4.set_xticks([0, 1, 2, 3, 4])
ax4.set_xticklabels(['$10^{-9}$', '$10^{-8}$', '$10^{-7}$', '$10^{-6}$', '$10^{-5}$'])
ax1.set_xticks([])
ax2.set_xticks([])

ax1.tick_params(labelsize = ticksize)
ax2.tick_params(labelsize = ticksize)
ax3.tick_params(labelsize = ticksize)
ax4.tick_params(labelsize = ticksize)

# Add colorbar
ax5 = fig.add_subplot(gs[1:, 3])
ax5.axis('off')
plt.text(1.2, -5.5, s=' $10^{-5}$ \n $10^{-6}$ \n $10^{-7}$ \n $10^{-8}$ \n $10^{-9}$', va='center', fontsize = ticksize, linespacing=10)
plt.text(3, -5.5, s='Variance of $[K^+]_i$ after 2000th pace', rotation = 90, fontsize = 35, va = 'center')
clb = plt.colorbar(image, cax=ax5, shrink = 0.3, ticks = [-8, -7, -6, -5, -4, -3])

# Add text labels
ax6 = fig.add_subplot(gs[:, 0])
ax6.axis('off')
ax6.text(x = -0.2, y = 0.65, s='TTP06', fontsize = 40, ha='center')
ax6.text(x = -0.2, y = 0.2, s='ORd CiPA', fontsize = 40, ha='center')

ax7 = fig.add_subplot(gs[0, :])
ax7.axis('off')
ax7.text(x = 0.35,y = 0, s='Algebraic voltage', fontsize = 40, ha='center')
ax7.text(x = 0.75, y = 0, s='Derivative voltage', fontsize = 40, ha='center')

# Save
plt.tight_layout(pad = 1.5)
plt.savefig('./Figures/Stability maps.png', dpi = 300)


# In[Plot the speed of simulation]
fig = plt.figure(constrained_layout=False, figsize = (24, 18))

gs = GridSpec(3, 4, figure=fig, width_ratios = [0.4, 1, 1, 0.1], height_ratios = [0.2, 1, 1],
              hspace = 0.1, wspace = 0.1)

# For the maps
ax1 = fig.add_subplot(gs[1, 1])
plt.imshow(durations_of_sim_an_TT06)
ax1.set_ylim([-0.5, 4.5])
plt.clim(0, 80)


ax2 = fig.add_subplot(gs[1, 2])
plt.imshow(durations_of_sim_der_TT06)
ax2.set_ylim([-0.5, 4.5])
plt.clim(0, 80)

ax3 = fig.add_subplot(gs[2, 1])
plt.imshow(durations_of_sim_an_ORd)
ax3.set_ylim([-0.5, 4.5])
plt.clim(0, 80)

ax4 = fig.add_subplot(gs[2, 2])
image = plt.imshow(durations_of_sim_der_ORd)
ax4.set_ylim([-0.5, 4.5])
plt.clim(0, 80)

# Label axes
ticksize = 22
labelsize = 28
ax1.set_ylabel('Relative solver tolerance', fontsize = labelsize)
ax3.set_ylabel('Relative solver tolerance', fontsize = labelsize)
ax1.set_yticks([0, 1, 2, 3, 4])
ax1.set_yticklabels(['$10^{-9}$', '$10^{-8}$', '$10^{-7}$', '$10^{-6}$', '$10^{-5}$'])
ax3.set_yticks([0, 1, 2, 3, 4])
ax3.set_yticklabels(['$10^{-9}$', '$10^{-8}$', '$10^{-7}$', '$10^{-6}$', '$10^{-5}$'])
ax2.set_yticks([])
ax4.set_yticks([])

ax3.set_xlabel('Absolute solver tolerance', fontsize = labelsize)
ax4.set_xlabel('Absolute solver tolerance', fontsize = labelsize)
ax3.set_xticks([0, 1, 2, 3, 4])
ax3.set_xticklabels(['$10^{-9}$', '$10^{-8}$', '$10^{-7}$', '$10^{-6}$', '$10^{-5}$'])
ax4.set_xticks([0, 1, 2, 3, 4])
ax4.set_xticklabels(['$10^{-9}$', '$10^{-8}$', '$10^{-7}$', '$10^{-6}$', '$10^{-5}$'])
ax1.set_xticks([])
ax2.set_xticks([])

ax1.tick_params(labelsize = ticksize)
ax2.tick_params(labelsize = ticksize)
ax3.tick_params(labelsize = ticksize)
ax4.tick_params(labelsize = ticksize)

# Add colorbar
ax5 = fig.add_subplot(gs[1:, 3])
ax5.axis('off')
plt.text(1.2, 40, s=' 80 \n 60 \n 40 \n 20 \n 0 ', va='center', fontsize = ticksize, linespacing=13)
plt.text(3, 40, s='Time of simulation (s)', rotation = 90, fontsize = 35, va = 'center')
clb = plt.colorbar(image, cax=ax5, shrink = 0.3, ticks = [-8, -7, -6, -5, -4, -3])

# Add text labels
ax6 = fig.add_subplot(gs[:, 0])
ax6.axis('off')
ax6.text(x = -0.2, y = 0.65, s='TTP06', fontsize = 40, ha='center')
ax6.text(x = -0.2, y = 0.2, s='ORd CiPA', fontsize = 40, ha='center')

ax7 = fig.add_subplot(gs[0, :])
ax7.axis('off')
ax7.text(x = 0.35,y = 0, s='Algebraic voltage', fontsize = 40, ha='center')
ax7.text(x = 0.75, y = 0, s='Derivative voltage', fontsize = 40, ha='center')

# Save
plt.tight_layout(pad = 1.5)
plt.savefig('./Figures/Speed maps.png', dpi = 300)

