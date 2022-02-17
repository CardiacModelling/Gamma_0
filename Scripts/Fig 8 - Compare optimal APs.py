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
os.chdir('C:/Users/barraly/Downloads/initial_conditions_study-master/initial_conditions_study-master')


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
    s.pre(2000000)
    out = s.run(1000, log_interval = 1)
    print('Potassium at steady-state: ' + str(np.round(out['intracellular_ions.ki'][-1], decimals = 2)))
    
    return np.array(out['membrane.V'])


# In[Reuse the fitting instructions]
# Define the time points on which to read the voltage
time_points = np.linspace(0, 999, 1000)

# Define the fitted parameters and initial point
parameters_to_fit = ['ical.rescale', 'ikr.rescale', 'IKs.rescale', 'INa.rescale', 'INaL.rescale']
true_values = np.array([1, 1, 1, 1, 1])


# In[Compute the fitted data]
# Set the parameters values
for p, label in enumerate(parameters_to_fit):
    s.set_constant(label, true_values[p])
# Run the model with the published original initial conditions and the Gamma_0 value associated with it
Gamma_0_for_fitting = -7.80116
data_to_fit = compute(Gamma_0_for_fitting)

# For validation with 50% IKr inhibition
s.set_constant(parameters_to_fit[1], 0.5 * true_values[1])
validation_data = compute(Gamma_0_for_fitting)


# In[Report the results from the fitting with Ohara initial conditions]
default_state = Ohara_init_conds.copy()
found_parameters = [1.000, 1.000, 1.000, 1.000, 1.000]
for p, label in enumerate(parameters_to_fit):
    s.set_constant(label, found_parameters[p])
fitting_with_Ohara_ICs = compute(Gamma_0_for_fitting)


print('Gamma_0 for fitting : ' + str(Gamma_0_for_fitting))
APD90 = sabs_pkpd.cardiac.compute_APD(AP = fitting_with_Ohara_ICs, time_points= time_points, upstroke_time = 50)
print('APD_90 at baseline : ' + str(APD90))

# Predict IKr block AP
s.set_constant(parameters_to_fit[1], 0.5 * found_parameters[1])
fitting_with_Ohara_ICs_Kr_blocked = compute(Gamma_0_for_fitting)
APD90 = sabs_pkpd.cardiac.compute_APD(AP = fitting_with_Ohara_ICs_Kr_blocked, time_points= time_points, upstroke_time = 50)
print('APD_90 after 50% IKr block : ' + str(APD90))


# In[Report the results from the fitting with TT06 initial conditions]
default_state[1] = 10.134
default_state[2] = 10.134
default_state[3] = 135.369
default_state[4] = 135.369
default_state[5] = 1.058e-04
default_state[6] = 2.142e-04
default_state[7] = 3.556
default_state[8] = 3.556
initial_voltage = -84.936

Gamma_0_for_fitting = G0_calc(Nai = default_state[1],
                              Nass = default_state[2],
                              Ki = default_state[3],
                              Kss = default_state[4],
                              Cai = default_state[5],
                              Cass = default_state[6],
                              Cansr = default_state[7],
                              Cajsr = default_state[8],
                              V=initial_voltage)

found_parameters = [0.8278844, 1.13276793, 0.74292672, 1.08754243, 1.4459109]
for p, label in enumerate(parameters_to_fit):
    s.set_constant(label, found_parameters[p])
fitting_with_TT06_ICs = compute(Gamma_0_for_fitting)

print('Gamma_0 for fitting : ' + str(Gamma_0_for_fitting))
APD90 = sabs_pkpd.cardiac.compute_APD(AP = fitting_with_TT06_ICs, time_points= time_points, upstroke_time = 50)
print('APD_90 at baseline : ' + str(APD90))

# Predict IKr block AP
s.set_constant(parameters_to_fit[1], 0.5 * found_parameters[1])
fitting_with_TT06_ICs_Kr_blocked = compute(Gamma_0_for_fitting)
APD90 = sabs_pkpd.cardiac.compute_APD(AP = fitting_with_TT06_ICs_Kr_blocked, time_points= time_points, upstroke_time = 50)
print('APD_90 after 50% IKr block : ' + str(APD90))


# In[Report the results from the fitting with TT06 initial conditions and Gamma_0]
default_state = Ohara_init_conds.copy()
found_parameters = [0.9999947, 0.99999936, 0.99995396, 0.999993485, 0.9999772, -7.801077]
for p, label in enumerate(parameters_to_fit):
    s.set_constant(label, found_parameters[p])
fitting_with_TT06_ICs_gamma_0 = compute(found_parameters[-1])


APD90 = sabs_pkpd.cardiac.compute_APD(AP = fitting_with_TT06_ICs_gamma_0, time_points= time_points, upstroke_time = 50)
print('APD_90 at baseline : ' + str(APD90))

# Predict IKr block AP
s.set_constant(parameters_to_fit[1], 0.5 * found_parameters[1])
fitting_with_TT06_ICs_gamma_0_Kr_blocked = compute(found_parameters[-1])
APD90 = sabs_pkpd.cardiac.compute_APD(AP = fitting_with_TT06_ICs_gamma_0_Kr_blocked, time_points= time_points, upstroke_time = 50)
print('APD_90 after 50% IKr block : ' + str(APD90))


# In[Plot the comparison]
def place_caption_label(ax, label, loc='upper left', fontsize=35):
    from matplotlib.offsetbox import AnchoredText
    at = AnchoredText(label, loc=loc, prop=dict(size=fontsize), frameon=True, borderpad = 0)
    ax.add_artist(at)    
    return None

from mpl_toolkits.axes_grid1.inset_locator import inset_axes
fig, ax = plt.subplots(1, 2, figsize=[15, 7])

size_ticks = 22
size_labels = 25
# Shift x-axis so that AP starts at 0 ms (in the simulations, the stimulus fires at t=50 ms)
x = np.linspace(0, 599, 600)

# Plot the fitted APs
ax[0].plot(x, data_to_fit[50:650], label = 'Data to fit', color = 'k', linewidth = 5)
ax[0].plot(x, fitting_with_Ohara_ICs[50:650], label = 'Fitting #1', linestyle = '--', linewidth = 3)
ax[0].plot(x, fitting_with_TT06_ICs[50:650], label = 'Fitting #2', linestyle = '--', linewidth = 3)
ax[0].plot(x, fitting_with_TT06_ICs_gamma_0[50:650], label = 'Fitting #3', linestyle = '--', linewidth = 3)

ax[0].legend(fontsize = 25)
ax[0].set_xlabel('Time (ms)', fontsize = size_labels)
ax[0].set_ylabel('Voltage (mV)', fontsize = size_labels)
ax[0].tick_params(axis = 'both', labelsize = size_ticks)
place_caption_label(ax[0], 'A', 'lower right')

# Add an inset to zoom into the short pacing periods
axins1 = ax[0].inset_axes(bounds = [0.7, 0.2, 0.3, 0.3])
x_inset = np.linspace(270, 299, 30)

axins1.plot(x_inset, data_to_fit[320:350], color = 'k', linewidth = 5)
axins1.plot(x_inset, fitting_with_Ohara_ICs[320:350], linestyle = '--', linewidth = 3)
axins1.plot(x_inset, fitting_with_TT06_ICs[320:350], linestyle = '--', linewidth = 3)
axins1.plot(x_inset, fitting_with_TT06_ICs_gamma_0[320:350], linestyle = '--', linewidth = 3)

# set up the inset ticks
axins1.set_xticks([270, 285, 300])
axins1.tick_params(axis = 'both', labelsize = 15)


# Plot the predicted APs with Kr block
ax[1].plot(x, validation_data[50:650], label = 'Validation data', linestyle = '-', linewidth = 5, color = 'k')
ax[1].plot(x, fitting_with_Ohara_ICs_Kr_blocked[50:650], label = 'Prediction #1', linestyle = '-', linewidth = 3)
ax[1].plot(x, fitting_with_TT06_ICs_Kr_blocked[50:650], label = 'Prediction #2', linestyle = '-', linewidth = 3)
ax[1].plot(x, fitting_with_TT06_ICs_gamma_0_Kr_blocked[50:650], label = 'Prediction #3', linestyle = '-', linewidth = 3)

ax[1].legend(fontsize = 25, loc = 'upper right')
ax[1].set_xlabel('Time (ms)', fontsize = size_labels)
ax[1].set_ylabel('Voltage (mV)', fontsize = size_labels)
ax[1].tick_params(axis = 'both', labelsize = size_ticks)
place_caption_label(ax[1], 'B', 'lower right')


# Add an inset to zoom into the short pacing periods
axins2 = ax[1].inset_axes(bounds = [0.7, 0.2, 0.3, 0.3])
x_inset = np.linspace(375, 424, 50)

axins2.plot(x_inset, validation_data[425:475], linestyle = '-', linewidth = 5, color = 'k')
axins2.plot(x_inset, fitting_with_Ohara_ICs_Kr_blocked[425:475], linestyle = '-', linewidth = 3)
axins2.plot(x_inset, fitting_with_TT06_ICs_Kr_blocked[425:475], linestyle = '-', linewidth = 3)
axins2.plot(x_inset, fitting_with_TT06_ICs_gamma_0_Kr_blocked[425:475], linestyle = '-', linewidth = 3)

# set up the inset ticks
axins2.set_xticks([375, 400, 425])
axins2.tick_params(axis = 'both', labelsize = 15)

# Save
plt.tight_layout()
plt.savefig('./Figures/Comparison of optimal APs.png', dpi = 300)



