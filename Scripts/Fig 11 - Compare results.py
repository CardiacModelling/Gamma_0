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
os.chdir('The/location/of/the/root/folder/of/this/repo')


# In[Load the model]
filename = './Models/Ohara CiPA - analytical voltage.mmt'
s = sabs_pkpd.load_model.load_simulation_from_mmt(filename)
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

def compute(Gamma_0, duration = 1000):
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
    out = s.run(duration, log_interval = 1)
    print('Potassium at steady-state: ' + str(np.round(out['intracellular_ions.ki'][-1], decimals = 2)))
    
    return np.array(out['membrane.V'])



# In[Get EADs with changed sodium initial conditions]
# Add ikr block to have EAD
s.set_constant('drug.ikr_rescale', 0.05)


# Run the model with the published original initial conditions and the Gamma_0 value associated with it
Gamma_0_for_EADs = -20

# Set the Nai to a high value
default_state[1] = 15
high_nai = compute(Gamma_0_for_EADs, duration = 2000)

# Redo the same with Nai with a low value
default_state[1] = 5
low_nai = compute(Gamma_0_for_EADs, duration = 2000)

# Original Nai intial value
default_state[1] =7.26800449799999981
original_nai = compute(Gamma_0_for_EADs, duration = 2000)


# In[Compute the fitted data]

# Run the model with the published original initial conditions and the Gamma_0 
# value associated with it
Gamma_0_for_EADs = -20

# Set the Nai to a high value
default_state[1] = 15

# Set the true value of parameters
parameters_to_fit = ['ical.rescale', 'ikr.rescale', 'IKs.rescale', 'INa.rescale', 'INaL.rescale', 'membrane.c0']
found_parameters = np.array([0, 0, 0, 0, 0, Gamma_0_for_EADs])
for p, label in enumerate(parameters_to_fit[:-1]):
    s.set_constant(label, np.exp(found_parameters[p]))
    
# Add ikr block to have EAD
s.set_constant('drug.ikr_rescale', 1)
data = compute(Gamma_0_for_EADs)
print('\nGamma_0 : ' + str(Gamma_0_for_EADs))
APD90 = sabs_pkpd.cardiac.compute_APD(AP = data, time_points= np.linspace(0, 999, 1000), upstroke_time = 50)
print('APD_90 at baseline : ' + str(APD90))

# Add ikr block to have EAD
s.set_constant('drug.ikr_rescale', 0.05)
data_block = compute(Gamma_0_for_EADs)
APD90 = sabs_pkpd.cardiac.compute_APD(AP = data_block, time_points= np.linspace(0, 999, 1000), upstroke_time = 50)
print('APD_90 at 95% IKr block : ' + str(APD90))
print('Diastolic K+ (95% IKr block) : ' + str(s.state()[3]))


# In[Reuse the fitting instructions]
# Define the time points on which to read the voltage
time_points = np.linspace(0, 999, 1000)

# Define the fitted parameters and initial point
parameters_to_fit = ['ical.rescale', 'ikr.rescale', 'IKs.rescale', 'INa.rescale', 'INaL.rescale', 'membrane.c0']
found_parameters = np.array([-1.46892781275010265e-01,-6.96374222510645069e-02,2.33739218682271238e-01,-6.60565445013814312e-02,4.53892759398383194e-01,-1.97465356401935992e+01])
for p, label in enumerate(parameters_to_fit[:-1]):
    s.set_constant(label, np.exp(found_parameters[p]))
    print(label + ' : ' + str(np.round(np.exp(found_parameters[p]), decimals = 3)))
# Reset the intial Nai to its original value used during the fitting
default_state[1] =7.26800449799999981

# Compute APD90 at baseline
s.set_constant('drug.ikr_rescale', 1)
fitted = compute(found_parameters[-1])
print('\nGamma_0 fitted : ' + str(found_parameters[-1]))
APD90 = sabs_pkpd.cardiac.compute_APD(AP = fitted, time_points= time_points, upstroke_time = 50)
print('APD_90 at baseline : ' + str(APD90))

# Compute APD90 at 95% IKr block
s.set_constant('drug.ikr_rescale', 0.05)
fitted_block = compute(found_parameters[-1])
APD90 = sabs_pkpd.cardiac.compute_APD(AP = fitted_block, time_points= time_points, upstroke_time = 50)
print('APD_90 at 95% IKr block : ' + str(APD90))
print('Diastolic K+ (95% IKr block) : ' + str(s.state()[3]))


# In[Plot the comparison]
def place_caption_label(ax, label, loc='upper left', fontsize=35):
    from matplotlib.offsetbox import AnchoredText
    at = AnchoredText(label, loc=loc, prop=dict(size=fontsize), frameon=True, borderpad = 0)
    ax.add_artist(at)    
    return None

fig, ax = plt.subplots(1, 2, figsize = (16, 8))

# Plot the training quality with Kr block
ax[0].plot(np.linspace(0, 999, 1000), data_block, label = 'True model', color = 'k', linewidth = 5)
ax[0].plot(np.linspace(0, 999, 1000), fitted_block, label = 'Fitted model', linestyle = '-', linewidth = 5)

# Labels
ax[0].legend(fontsize = 20, loc= 'upper right')
ax[0].set_xlabel('Time (ms)', fontsize = 28)
ax[0].set_ylabel('Voltage (mV)', fontsize = 28)
ax[0].tick_params(labelsize= 25)
place_caption_label(ax[0], 'A', loc='lower right', fontsize=35)
ax[0].text(x=800, y = -10, s = 'Model\ntraining', ha='center', fontsize = 40)

# Plot the fitted APs
ax[1].plot(np.linspace(0, 499, 500), data[:500], label = 'True model, baseline AP', color = 'k', linewidth = 5)
ax[1].plot(np.linspace(0, 499, 500), fitted[:500], label = 'Fitted model, baseline AP', linestyle = '-', linewidth = 5)

ax[1].set_xlabel('Time (ms)', fontsize = 28)
ax[1].set_ylabel('Voltage (mV)', fontsize = 28)
ax[1].set_xticks([0, 200, 400])
ax[1].tick_params(labelsize = 25)
place_caption_label(ax[1], 'B', loc='lower right', fontsize=35)
ax[1].text(x=400, y = -10, s = 'Validation\nprediction', ha='center', fontsize = 40)

# Save
plt.tight_layout()
plt.savefig('./Figures/Optimisation fitted with EAD bifurcation.png', dpi = 300)