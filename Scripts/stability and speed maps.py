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


# In[Load the models for TTP06]
# TTP06 with analytical voltage expression
analytic = sabs_pkpd.load_model.load_simulation_from_mmt('C:/Users/barraly/Documents/PhD/MMT models/tentusscher_2006 - analytical voltage.mmt')

# TTP06 with derivative voltage expression
derivative = sabs_pkpd.load_model.load_simulation_from_mmt('C:/Users/barraly/Documents/PhD/MMT models/tentusscher_2006.mmt')

readout = 'potassium.K_i'


# In[Load the models for ORd]
# ORd with analytical voltage expression
analytic = sabs_pkpd.load_model.load_simulation_from_mmt('C:/Users/barraly/Documents/PhD/MMT models/Ohara CiPA - analytical voltage.mmt')

# ORd with derivative voltage expression
derivative = sabs_pkpd.load_model.load_simulation_from_mmt('C:/Users/barraly/Documents/PhD/MMT models/ohara_rudy_cipa_v1_2017.mmt')

readout = 'intracellular_ions.ki'


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
    plt.xlabel('# Paces', Fontsize = 17)
    plt.ylabel('$K_i$ (mM)', Fontsize = 17)
    plt.title('Solver tolerance : ' + str(tol), Fontsize = 24)
    
    plt.subplot(1, 2, 2)
    plt.plot(np.linspace(2000, 3000, 1000), an[readout][2000:], label = 'Analytical')
    plt.plot(np.linspace(2000, 3000, 1000), der[readout][2000:], label = 'Derivative')
    plt.legend(fontsize = 15)
    plt.xlabel('# Paces', Fontsize = 17)
    plt.ylabel('$K_i$ (mM)', Fontsize = 17)
    plt.title('Solver tolerance : ' + str(tol), Fontsize = 24)
    
    plt.tight_layout()
    

# In[Try different solver tolerances]
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
        
               
# In[ Plot the convergence to limit cycle and stability after the 2000th pace]

plt.figure(figsize = (12, 9))
plt.subplot(1, 2, 1)
plt.title('Analytical voltage', Fontsize = 20)
plt.imshow(np.log10(stability_analytic))
plt.xlabel('Absolute solver tolerance', Fontsize = 17)
plt.xticks(ticks = [0, 1, 2, 3, 4], labels = ['1e-09', '1e-08', '1e-07', '1e-06', '1e-05'])
plt.ylim([-0.5, 4.5])
plt.ylabel('Relative solver tolerance', Fontsize = 17)
plt.yticks(ticks = [0, 1, 2, 3, 4], labels = ['1e-09', '1e-08', '1e-07', '1e-06', '1e-05'])
plt.clim(-8.5, -2.5)
plt.colorbar(shrink = 0.5)
plt.text(5.7, -0.5, s='Log of mean distance to 2000th pace', rotation = 90, Fontsize = 15)

plt.subplot(1, 2, 2)
plt.title('Derivative voltage', Fontsize = 20)
plt.imshow(np.log10(stability_derivative))
plt.xlabel('Absolute solver tolerance', Fontsize = 17)
plt.xticks(ticks = [0, 1, 2, 3, 4], labels = ['1e-09', '1e-08', '1e-07', '1e-06', '1e-05'])
plt.ylim([-0.5, 4.5])
plt.ylabel('Relative solver tolerance', Fontsize = 17)
plt.yticks(ticks = [0, 1, 2, 3, 4], labels = ['1e-09', '1e-08', '1e-07', '1e-06', '1e-05'])
plt.clim(-8.5, -2.5)

plt.colorbar(shrink = 0.5)
plt.text(5.7, -0.5, s='Log of mean distance to 2000th pace', rotation = 90, Fontsize = 15)

plt.tight_layout(pad = 2.5)


# In[Plot the speed of simulation]
plt.figure(figsize = (12, 9))
plt.subplot(1, 2, 1)
plt.title('Analytical voltage', Fontsize = 20)
plt.imshow(durations_of_sim_an)
plt.xlabel('Absolute solver tolerance', Fontsize = 17)
plt.xticks(ticks = [0, 1, 2, 3, 4], labels = ['1e-09', '1e-08', '1e-07', '1e-06', '1e-05'])
plt.ylim([-0.5, 4.5])
plt.ylabel('Relative solver tolerance', Fontsize = 17)
plt.yticks(ticks = [0, 1, 2, 3, 4], labels = ['1e-09', '1e-08', '1e-07', '1e-06', '1e-05'])
plt.clim(0, 30)
plt.colorbar(shrink = 0.5)
plt.text(5.7, 0, s='Duration of simulation (in s)', rotation = 90, Fontsize = 15)

plt.subplot(1, 2, 2)
plt.title('Derivative voltage', Fontsize = 20)
plt.imshow(durations_of_sim_der)
plt.xlabel('Absolute solver tolerance', Fontsize = 17)
plt.xticks(ticks = [0, 1, 2, 3, 4], labels = ['1e-09', '1e-08', '1e-07', '1e-06', '1e-05'])
plt.ylim([-0.5, 4.5])
plt.ylabel('Relative solver tolerance', Fontsize = 17)
plt.yticks(ticks = [0, 1, 2, 3, 4], labels = ['1e-09', '1e-08', '1e-07', '1e-06', '1e-05'])
plt.clim(0, 30)

plt.colorbar(shrink = 0.5)
plt.text(5.7, 0, s='Duration of simulation (in s)', rotation = 90, Fontsize = 15)

plt.tight_layout(pad = 2.5)

