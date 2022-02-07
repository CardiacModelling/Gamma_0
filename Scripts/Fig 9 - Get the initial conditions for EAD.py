# -*- coding: utf-8 -*-
"""
Created on Tue Mar  9 09:42:00 2021

@author: barraly
"""

import sabs_pkpd
import numpy as np
import matplotlib.pyplot as plt
import os

# Select the folder in which this repo is downloaded in the line below
# os.chdir('The/location/of/the/root/folder/of/this/repo')
os.chdir('C:/Users/barraly/Downloads/initial_conditions_study-master/initial_conditions_study-master')


# In[Load the model]
s = sabs_pkpd.load_model.load_simulation_from_mmt('./Models/Ohara CiPA - analytical voltage.mmt')
s.set_tolerance(1e-08, 1e-08)
default_state = s.state()

# Save the initial conditions published in OHara CiPA model
Ohara_init_conds = default_state.copy()


# In[Defin the needed functions]
# Define the functions to make sure there is consistency between the initial conditions
def G0_calc(Ki = 144.65559, Kss = 144.65556, Nai = 7.268, Nass = 7.26809,
            Cai = 8.6e-5, Cansr = 1.61957, Cajsr = 1.571234014, Cass = 8.49e-5,
            V=-88, extraK = 5.4, extraNa = 140, extraCa = 1.8):
    tot_cai = Cai * (1 + 0.05 / (Cai + 0.00238) + 0.07/(Cai + 0.0005))
    tot_cass = Cass * (1 + 0.047 / (Cass + 0.00087) + 1.124/(Cass + 0.0087))
    tot_cajsr = Cajsr * (1 + 10 / (Cajsr + 0.8))
    return V / (96485 * 2.583592e-05) * 0.0001533576 - (Ki + Kss * 0.029411764705882353 + Nai + Nass * 0.029411764705882353 + 2*(tot_cai + tot_cass * 0.029411764705882353 + Cansr * 0.08117647059 + tot_cajsr * 0.007059) - extraK - extraNa - 2 * extraCa)

def Ki_calc(G0, Nai = 7.268, Nass = 7.26809, Cai = 8.6e-5, Cansr = 1.61957, Cajsr = 1.571234014, Cass = 8.49e-5, V=-88, extraK = 5.4, extraNa = 140, extraCa = 1.8):
    tot_cai = Cai * (1 + 0.05 / (Cai + 0.00238) + 0.07/(Cai + 0.0005))
    tot_cass = Cass * (1 + 0.047 / (Cass + 0.00087) + 1.124/(Cass + 0.0087))
    tot_cajsr = Cajsr * (1 + 10 / (Cajsr + 0.8))
    return (V / (96485 * 2.583592e-05) * 0.0001533576 + extraK + extraNa + 2 * extraCa - G0 - Nai - Nass * 0.029411764705882353 - 2*(tot_cai + tot_cass * 0.029411764705882353 + Cansr * 0.08117647059 + tot_cajsr * 0.007059)) / 1.029411764705882353 

def compute(Gamma_0):
    # Reinitialise the myokit.Simulation
    s.reset()
    
    # Set the initial conditions for Ki and Kss so that the initial conditions match with the value of Gamma_0
    initial_state = default_state.copy()
    initial_K = Ki_calc(Gamma_0,
                       Nai = default_state[1],
                       Nass = default_state[2],
                       Cai = default_state[5],
                       Cansr = default_state[7],
                       Cajsr = default_state[8],
                       Cass = default_state[6])
    initial_state[3] = initial_K
    initial_state[4] = initial_K
    s.set_state(initial_state)
    
    # Set the value of Gamma_0 in the myokit.Simulation
    s.set_constant('membrane.c0', Gamma_0)
    
    # Record the action potential at the limit cycle
    s.pre(1500000)
    out = s.run(2000, log_interval = 1)
    print('Potassium at steady-state: ' + str(np.round(out['intracellular_ions.ki'][-1], decimals = 2)))
    
    return np.array(out['membrane.V'])


# In[Get EADs with changed sodium initial conditions]
# Add ikr block to have EAD
s.set_constant('drug.ikr_rescale', 0.05)


# Run the model with the published original initial conditions and the Gamma_0 value associated with it
Gamma_0_for_EADs = -20

# Set the Nai to a high value
default_state[1] = 15
high_nai = compute(Gamma_0_for_EADs)

# Redo the same with Nai with a low value
default_state[1] = 5
low_nai = compute(Gamma_0_for_EADs)

# Original Nai intial value
default_state[1] =7.26800449799999981
original_nai = compute(Gamma_0_for_EADs)


# In[Plot results]
plt.figure(figsize = (12, 8))
plt.plot(high_nai, label = 'Initial $[Na^+]_i = 15 mM$', linewidth = 5, color = 'k')
#plt.plot(low_nai, label = 'Low initial Nai', linewidth = 3)
plt.plot(original_nai, label = 'Initial $[Na^+]_i = 7.3 mM$', linewidth = 5)
plt.legend(fontsize = 26)
plt.xlabel('Time (ms)', fontsize = 29)
plt.ylabel('Voltage (mV)', fontsize = 29)
plt.xticks([0, 500, 1000, 1500, 2000], [0, 500, 1000, 1500, 2000], fontsize = 26)
plt.yticks(fontsize = 26)
#plt.title('Steady-state APs with 95% block of $I_{Kr}$\n with $\Gamma_0 = -20$ and varying initial $Na_i$ values', fontsize = 25)

# Save
plt.tight_layout()
plt.savefig('./Figures/bifurcation APs.png', dpi = 300)