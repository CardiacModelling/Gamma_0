# -*- coding: utf-8 -*-
"""
Created on Tue May 12 14:25:17 2020

@author: yann-stanislas.barral@roche.com
"""

import numpy as np
import sabs_pkpd
import matplotlib.pyplot as plt
import matplotlib.axes
from mpl_toolkits.mplot3d import Axes3D


# In[Load the models]
s = sabs_pkpd.load_model.load_simulation_from_mmt('./Ohara CiPA - analytical voltage - paper.mmt')
s1 = sabs_pkpd.load_model.load_simulation_from_mmt('./tentusscher_2006 - analytical voltage.mmt')
default = s.state()
default1 = s1.state()


# In[For ORd CiPA]
# Define the functions to compute G0 and Ki from the other variables
def gamma0_calc(Ki = 144.65559, Kss = 144.65556, Nai = 7.268, Nass = 7.26809,
            Cai = 8.6e-5, Cansr = 1.61957, Cajsr = 1.571234014, Cass = 8.49e-5,
            V=-88, extraK = 5.4, extraNa = 140, extraCa = 1.8):
    tot_cai = Cai * (1 + 0.05 / (Cai + 0.00238) + 0.07/(Cai + 0.0005))
    tot_cass = Cass * (1 + 0.047 / (Cass + 0.00087) + 1.124/(Cass + 0.0087))
    tot_cajsr = Cajsr * (1 + 10 / (Cajsr + 0.8))
    return -V / (96485 * 2.583592e-05) * 0.0001533576 + Ki + Kss * 0.029411764705882353 + Nai + Nass * 0.029411764705882353 + 2*(tot_cai + tot_cass * 0.029411764705882353 + Cansr * 0.08117647059 + tot_cajsr * 0.007059) - extraK - extraNa - extraCa

def Ki_calc(G0, Nai = 7.268, Nass = 7.26809, Cai = 8.6e-5, Cansr = 1.61957, Cajsr = 1.571234014, Cass = 8.49e-5, V=-88, extraK = 5.4, extraNa = 140, extraCa = 1.8):
    tot_cai = Cai * (1 + 0.05 / (Cai + 0.00238) + 0.07/(Cai + 0.0005))
    tot_cass = Cass * (1 + 0.047 / (Cass + 0.00087) + 1.124/(Cass + 0.0087))
    tot_cajsr = Cajsr * (1 + 10 / (Cajsr + 0.8))
    return (V / (96485 * 2.583592e-05) * 0.0001533576 + extraK + extraNa + extraCa + G0 - Nai - Nass * 0.029411764705882353 - 2*(tot_cai + tot_cass * 0.029411764705882353 + Cansr * 0.08117647059 + tot_cajsr * 0.007059)) / 1.029411764705882353 

# Define the function to run the model to steady-state from 
# a given starting point
def run_sim(concs, fixed_G0):
    # Retrieve the initial values for the state variables
    Nai = concs[0]
    Ki = concs[1]
    Cai = concs[2]
    Cass = concs[3]
    Cansr = concs[4]
    Cajsr = concs[5]
    
    # Distinguish cases where G0 is fixed throughout the simulations
    if fixed_G0 == True:
        # Compute Ki to stay on the same G0 hyperplane
        Ki = Ki_calc(G0, Nai, Nai, Cai, Cansr, Cajsr, Cass)
    
    else:
        G0_sim = gamma0_calc(Ki, Ki, Nai, Nai, Cai, Cansr, Cajsr, Cass)
        s.set_constant('membrane.c0', G0_sim)
        
    # rescale of initial conditions
    init_conds = default[:]
    init_conds[1] = Nai
    init_conds[2] = Nai
    init_conds[3] = Ki
    init_conds[4] = Ki
    init_conds[5] = Cai
    init_conds[6] = Cass
    init_conds[7] = Cansr
    init_conds[8] = Cajsr
    
    # Generate random initial conditions for the other state variables
    # Note that the Markovian states of IKr model must sum up to 1
    init_conds [9:] = np.random.uniform(0, 1, size = (len(init_conds[9:])))
    init_conds[33] = 1 - init_conds[34] - init_conds[35] - init_conds[36] - init_conds[37] - init_conds[38]

    # Set the bound IKr variables to 0
    init_conds[39] = 0
    init_conds[41] = 0
    init_conds[40] = 0
            
    # Set the initial conditions to the default state and reset the rest
    s.set_default_state(init_conds)
    s.reset()
    
    # Run to steady-state wihtout pacing
    # Run first for 500 ms to move past the first peak and enable smooth plotting of convergence
    adapt = s.run(500)
    a = s.run(5000000, log=['intracellular_ions.ki', 'intracellular_ions.nai', 'membrane.V'])
    
    return a

# Set solver tolerance to fine
s.set_tolerance(abs_tol=1e-08, rel_tol = 1e-06)

# Set no pacing
s.set_constant('stimulus_protocol.i_Stim_Amplitude', 0)

# Set G0 to the value that is supposed to be constant
G0 = gamma0_calc(Ki = default[3], Kss = default[4], Nai = default[1], Nass = default[2],
             Cai = default[5], Cansr = default[7], Cajsr = default[8], Cass = default[6])
s.set_constant('membrane.c0', G0)

# Prepare the sampling points of the G0 hyperplane
nb_points = 30
points = np.zeros((nb_points, 6))
# Nais to be sampled from the physiological range of concentrations
points[:, 0] = np.random.uniform(4, 16, size = nb_points)
# Cais
points[:, 2] = np.random.uniform(0.5, 1.5, size = nb_points) * default[5]
#Cass
points[:, 3] = np.random.uniform(0.5, 1.5, size = nb_points) * default[6]
#Cansr
points[:, 4] = np.random.uniform(0.5, 1.5, size = nb_points) * default[7]
#Cajsr
points[:, 5] = np.random.uniform(0.5, 1.5, size = nb_points) * default[8]

# Run the simulations
outputs = []
for i in range(nb_points):
    outputs.append(run_sim(points[i, :], fixed_G0 = True))
    print(i)
    
# Redo the same with variation of all of the ICs independently, allowing variations of C_0

# Prepare the sampling points of the G0 hyperplane
points_1 = np.zeros((nb_points, 6))
# Nais to be sampled from the physiological range of concentrations
points_1[:, 0] = np.random.uniform(4, 16, size = nb_points)
# Kis to be sampled from the physiological range of concentrations
points_1[:, 1] = np.random.uniform(120, 153, size = nb_points)
# Cais
points_1[:, 2] = np.random.uniform(0.5, 1.5, size = nb_points) * default[5]
#Cass
points_1[:, 3] = np.random.uniform(0.5, 1.5, size = nb_points) * default[6]
#Cansr
points_1[:, 4] = np.random.uniform(0.5, 1.5, size = nb_points) * default[7]
#Cajsr
points_1[:, 5] = np.random.uniform(0.5, 1.5, size = nb_points) * default[8]

# Run the simulations
outputs_var_G0 = []
for i in range(nb_points):
    outputs_var_G0.append(run_sim(points_1[i, :], fixed_G0 = False))
    print(i)



# In[Plot the results for ORd CiPA]
Nais = []
Kis = []
V = []

for i in range(30):
    Nais.append(outputs[i]['intracellular_ions.nai'])
    Kis.append(outputs[i]['intracellular_ions.ki'])
    V.append(outputs[i]['membrane.V'])

fig = plt.figure(figsize = (13, 13))
ax = fig.add_subplot(111, projection='3d')

# Plot the evolution of the state variables versus time
for i in range(30):
    ax.plot(Nais[i], Kis[i], V[i])

# Plot the starting points
for i in range(30):
    ax.scatter(Nais[i][0], Kis[i][0], V[i][0], marker = '+', LineWidth = 20)
    
# Plot the state varaibles at quiescent steady-state
for i in range(30):
    ax.scatter(Nais[i][-1], Kis[i][-1], V[i][-1], LineWidth  = 10)

ax.view_init(elev=20, azim=-160)
ax.set_xlabel('Nai (mM)', Fontsize = 30, labelpad = 15)
ax.set_ylabel('Ki (mM)', Fontsize = 30, labelpad = 15)
ax.set_zlabel('Voltage (mV)', Fontsize = 30, labelpad = 15)
ax.set_zlim([-89, -87.5])



Nais = []
Kis = []
V = []

for i in range(30):
    Nais.append(outputs_var_G0[i]['intracellular_ions.nai'])
    Kis.append(outputs_var_G0[i]['intracellular_ions.ki'])
    V.append(outputs_var_G0[i]['membrane.V'])

fig = plt.figure(figsize = (13, 13))
ax = fig.add_subplot(111, projection='3d')

# Plot the evolution of the state variables versus time
for i in range(30):
    ax.plot(Nais[i], Kis[i], V[i])
    
# Plot the starting points
for i in range(30):
    ax.scatter(Nais[i][0], Kis[i][0], V[i][0], marker = '+', LineWidth = 20)
    
# Plot the state varaibles at quiescent steady-state
for i in range(30):
    ax.scatter(Nais[i][-1], Kis[i][-1], V[i][-1], LineWidth  = 10)

ax.view_init(elev=20, azim=-160)
ax.set_xlabel('Nai (mM)', Fontsize = 30, labelpad = 15)
ax.set_ylabel('Ki (mM)', Fontsize = 30, labelpad = 15)
ax.set_zlabel('Voltage (mV)', Fontsize = 30, labelpad = 15)
ax.set_zlim([-90,  -83])


# In[For TT06]
# Define the functions to compute G0 and Ki from the other variables
def G0_calc(Ki, Nai, Cai = 0.00010578226860054304, Casr = 3.5556779165585235,
            Cass = 0.00021417287326980984, V=-85.3, extraK = 5.4, extraNa = 140,
            extraCa = 2):
    
    return -V / 96.4853415 * 0.000185 / 0.016404 + (Ki + Nai + 2*(Cai + Cai*0.2 / (Cai + 0.001)) +
                                                 2*1.094/16.404*(Casr + Casr * 10 / (Casr + 0.3)) +
                                                 2*5.468e-2/16.404*(Cass + Cass * 0.4 / (Cass + 0.00025)) - extraK - extraNa - extraCa)

def Ki_calc(G0, Nai, Cai, Casr, Cass, V = -85.3, extraK = 5.4, extraNa = 140,
            extraCa = 2):
    Cai_tot = Cai + Cai*0.2 / (Cai + 0.001)
    Casr_tot = Casr + Casr * 10 / (Casr + 0.3)
    Cass_tot = Cass + Cass * 0.4 / (Cass + 0.00025)
    Ki = extraK + extraNa + extraCa + G0 + V / 96.4853415 * 0.000185 / 0.016404 - (Nai + 2*Cai_tot  + 2*0.001094/0.016404*Casr_tot + 2*5.468e-5/0.016404*Cass_tot)
    return Ki

# Define the function to run the model to steady-state from 
# a given starting point
def run_sim(concs, fixed_G0):
    # Retrieve the initial values for the state variables
    Nai = concs[0]
    Ki = concs[1]
    Cai = concs[2]
    Casr = concs[3]
    Cass = concs[4]
    
    # rescale of initial conditions
    init_conds = default1[:]
    init_conds[0] = Cai
    init_conds[1] = Casr
    init_conds[2] = Cass
    init_conds[4] = Nai
    init_conds[5] = Ki
    
    # Distinguish cases where G0 is fixed or not
    if fixed_G0 == True:
        # Compute Ki to stay on the same G0 hyperplane
        Ki = Ki_calc(G0, Nai, Cai, Casr, Cass, V=-85.3)
        init_conds[5] = Ki
        # Set G0
        s1.set_constant('membrane.c0', G0)
        
    if fixed_G0 == False:
        # Compute Ki to stay on the same G0 hyperplane
        G0_sim = G0_calc(Ki, Nai, Cai, Casr, Cass, V=-85.3)
        s1.set_constant('membrane.c0', G0_sim)
    
    # Generate random initial conditions for the other state variables
    init_conds [6:] = np.random.uniform(0, 1, size = (len(init_conds[6:])))

    
    # Set the initial conditions to the default state and reset the rest
    s1.set_default_state(init_conds)
    s1.reset()
        
    # Set no pacing
    s1.set_constant('stimulus.amplitude', 0)

    # Run to steady-state wihtout pacing
    peak = s1.run(500)
    a = s1.run(5000000, log=['potassium.K_i', 'sodium.Na_i', 'membrane.V'])
    
    return a


# Set solver tolerance to fine
s1.set_tolerance(abs_tol=1e-08, rel_tol = 1e-06)

# Set G0 to the value that is supposed to be constant
G0 = G0_calc(Ki = default1[5], Nai = default1[4], Cai = default1[0], Casr = default1[1], Cass = default1[2], V = -85.3)

# Prepare the sampling points of the G0 hyperplane
nb_points = 30
points = np.zeros((nb_points, 5))
# Nais to be sampled from the physiological range of concentrations
points[:, 0] = np.random.uniform(4, 16, size = nb_points)
# Cai
points[:, 2] = np.random.uniform(0.5, 1.5, size = nb_points) * default1[0]
#Casr
points[:, 3] = np.random.uniform(0.5, 1.5, size = nb_points) * default1[1]
#Cass
points[:, 4] = np.random.uniform(0.5, 1.5, size = nb_points) * default1[2]

# Run the simulations
outputs1 = []
for i in range(nb_points):
    outputs1.append(run_sim(points[i, :], fixed_G0 = True))
    print(i)


# Redo the same with independent variations of all initial conditions, allowing variations of C_0    
# Prepare the sampling points of the G0 hyperplane
points1_1 = np.zeros((nb_points, 5))
# Nais to be sampled from the physiological range of concentrations
points1_1[:, 0] = np.random.uniform(4, 16, size = nb_points)
# Kis to be sampled from the physiological range of concentrations
points1_1[:, 1] = np.random.uniform(120, 145, size = nb_points)
# Cai
points1_1[:, 2] = np.random.uniform(0.5, 1.5, size = nb_points) * default1[0]
#Casr
points1_1[:, 3] = np.random.uniform(0.5, 1.5, size = nb_points) * default1[1]
#Cass
points1_1[:, 4] = np.random.uniform(0.5, 1.5, size = nb_points) * default1[2]

# Run the simulations
outputs1_var_G0 = []
for i in range(nb_points):
    outputs1_var_G0.append(run_sim(points1_1[i, :], fixed_G0 = False))
    print(i)


# In[Plot the results for TT06]
fig = plt.figure(figsize = (13, 13))
ax = fig.add_subplot(111, projection='3d')

# Plot the evolution of the state variables versus time
for i in range(30):
    ax.plot(outputs1[i]['sodium.Na_i'], outputs1[i]['potassium.K_i'], outputs1[i]['membrane.V'])

# Plot the starting points
for i in range(30):
    ax.scatter(outputs1[i]['sodium.Na_i'][0], outputs1[i]['potassium.K_i'][0], outputs1[i]['membrane.V'][0], marker = '+', LineWidth = 20)
   
# Plot the state varaibles at quiescent steady-state
for i in range(30):
    ax.scatter(outputs1[i]['sodium.Na_i'][-1], outputs1[i]['potassium.K_i'][-1], outputs1[i]['membrane.V'][-1], LineWidth  = 10)

ax.view_init(elev=20, azim=-160)
ax.set_xlabel('Nai (mM)', Fontsize = 30, labelpad = 15)
ax.set_ylabel('Ki (mM)', Fontsize = 30, labelpad = 15)
ax.set_zlabel('Voltage (mV)', Fontsize = 30, labelpad = 15)
ax.set_zlim([-87, -85])


fig = plt.figure(figsize = (13, 13))
ax = fig.add_subplot(111, projection='3d')

# Plot the evolution of the state variables versus time
for i in range(30):
    ax.plot(outputs1_var_G0[i]['sodium.Na_i'], outputs1_var_G0[i]['potassium.K_i'], outputs1_var_G0[i]['membrane.V'])

# Plot the starting points
for i in range(30):
    ax.scatter(outputs1_var_G0[i]['sodium.Na_i'][0], outputs1_var_G0[i]['potassium.K_i'][0], outputs1_var_G0[i]['membrane.V'][0], marker = '+', LineWidth = 20)
   
# Plot the state varaibles at quiescent steady-state
for i in range(30):
    ax.scatter(outputs1_var_G0[i]['sodium.Na_i'][-1], outputs1_var_G0[i]['potassium.K_i'][-1], outputs1_var_G0[i]['membrane.V'][-1], LineWidth  = 10)

ax.view_init(elev=20, azim=-160)
ax.set_xlabel('Nai (mM)', Fontsize = 30, labelpad = 15)
ax.set_ylabel('Ki (mM)', Fontsize = 30, labelpad = 15)
ax.set_zlabel('Voltage (mV)', Fontsize = 30, labelpad = 15)
ax.set_zlim([-90, -83])
