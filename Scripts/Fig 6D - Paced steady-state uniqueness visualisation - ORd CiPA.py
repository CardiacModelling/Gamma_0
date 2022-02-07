# -*- coding: utf-8 -*-
"""
Created on Tue Dec 8 16:04:17 2020

@author: yann-stanislas.barral@roche.com
"""

import numpy as np
import sabs_pkpd
import matplotlib.pyplot as plt
import matplotlib.axes
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.colors as colors
import matplotlib.cm as cm
import os

# Select the folder in which this repo is downloaded in the line below
os.chdir('The/location/of/the/root/folder/of/this/repo')



# In[Load the models]
s = sabs_pkpd.load_model.load_simulation_from_mmt('./Models/Ohara CiPA - analytical voltage.mmt')
default = s.state()

# Set solver tolerance to fine
s.set_tolerance(abs_tol=1e-08, rel_tol = 1e-08)


# In[Check that the model is close to steady-state]
# Reset the simulation
s.set_state(default)
s.set_time(0)

# Run the simulation to stead-state
convergence = s.run(2000000)
steady = s.run(1000)

# Plot the limit cycle
plt.figure(figsize = (7, 7))
plt.ticklabel_format(style='plain')
plt.plot(steady['environment.time'], steady['membrane.V'])
plt.xlabel('Time (ms)', fontsize  = 20)
plt.ylabel('Voltage (mV)', fontsize = 20)


# Plot the convergence towards equilibrium
plt.figure(figsize = (7, 7))
plt.plot(convergence['environment.time'], convergence['intracellular_ions.ki'])
plt.ylabel('$[K^+]_i$ (mM)', fontsize = 20)
plt.xlabel('Time (ms)', fontsize = 20)


# In[For ORd CiPA]
# Define the functions to compute G0 and Ki from the other variables
def G0_calc(Ki = 144.65559, Kss = 144.65556, Nai = 7.2680045, Nass = 7.26809,
            Cai = 8.6e-5, Cansr = 1.61957, Cajsr = 1.571234014, Cass = 8.49e-5,
            V=-88, Nao = 140, Ko = 5.4, Cao = 1.8):
    tot_cai = Cai * (1 + 0.05 / (Cai + 0.00238) + 0.07/(Cai + 0.0005))
    tot_cass = Cass * (1 + 0.047 / (Cass + 0.00087) + 1.124/(Cass + 0.0087))
    tot_cajsr = Cajsr * (1 + 10 / (Cajsr + 0.8))
    G0 = V / (96485 * 2.583592e-05) * 0.0001533576 - (Ki + Kss * 0.029411764705882353 + Nai + Nass * 0.029411764705882353 + 2*(tot_cai + tot_cass * 0.029411764705882353 + Cansr * 0.08117647059 + tot_cajsr * 0.007059) - Nao - 2*Cao - Ko)

def Ki_calc(G0, Nai = 7.268, Nass = 7.26809, Cai = 8.6e-5, Cansr = 1.61957, Cajsr = 1.571234014, Cass = 8.49e-5, V=-88, extraK = 5.4, extraNa = 140, extraCa = 1.8):
    tot_cai = Cai * (1 + 0.05 / (Cai + 0.00238) + 0.07/(Cai + 0.0005))
    tot_cass = Cass * (1 + 0.047 / (Cass + 0.00087) + 1.124/(Cass + 0.0087))
    tot_cajsr = Cajsr * (1 + 10 / (Cajsr + 0.8))
    return (V / (96485 * 2.583592e-05) * 0.0001533576 + extraK + extraNa + 2 * extraCa - G0 - Nai - Nass * 0.029411764705882353 - 2*(tot_cai + tot_cass * 0.029411764705882353 + Cansr * 0.08117647059 + tot_cajsr * 0.007059)) / 1.029411764705882353 

# Define the function to sample a G0 plane
def sample_hyperplane(G0, nb_points_per_axis = 15):
    nb_points = nb_points_per_axis**2
    points = np.zeros((nb_points, 8))
    
    for i in range(nb_points_per_axis):
        # Nai to be sampled from the physiological range of concentrations
        points[nb_points_per_axis * i:nb_points_per_axis * (i+1), 0] = np.linspace(4, 16, nb_points_per_axis)
        # Nass set to the same value initially as Nai
        points[nb_points_per_axis * i:nb_points_per_axis * (i+1), 1] = points[nb_points_per_axis * i:nb_points_per_axis * (i+1), 0]
        # Cai
        points[nb_points_per_axis * i:nb_points_per_axis * (i+1), 4] = np.linspace(0.5, 1.5, nb_points_per_axis) * default[5]
        #Cass
        points[nb_points_per_axis * i:nb_points_per_axis * (i+1), 5] = np.random.uniform(0.5, 1.5, size = nb_points_per_axis) * default[6]
        #Cansr
        points[nb_points_per_axis * i:nb_points_per_axis * (i+1), 6] = np.random.uniform(0.5, 1.5, size = nb_points_per_axis) * default[7]
        #Cajsr
        points[nb_points_per_axis * i:nb_points_per_axis * (i+1), 7] = np.random.uniform(0.5, 1.5, size = nb_points_per_axis) * default[8]
    
    # Compute Ki for each point now to stick to the same G0 value
    for sample in range(nb_points):
        points[sample, 2] = Ki_calc(G0,
              Nai = points[sample, 0],
              Nass = points[sample, 1],
              Cai = points[sample, 4],  
              Cass = points[sample, 5],                   
              Cansr = points[sample, 6],                   
              Cajsr = points[sample, 7],
              V=-88)
        # Set Kss to the same initial value as Ki
        points[sample, 3] = points[sample, 2]
        
    return points

# Define the function to run the model to steady-state from 
# a given starting point
def run_sim(concs):
    # Define the initial state with the intracellular concentrations provided
    init = default[:]
    init[1:9] = concs
    
    # Set non-concentrations state variables to random
    init[0] = np.random.uniform(0.5, 1.5) * init[0]
    init[9:] = np.random.uniform(0.0000001, 0.999999, size = 38)
    
    # Set IKr state variables to the value matching with a sum of states equal to 1
    init[33] = np.random.uniform(0.0000001, 0.999999)
    for i in range(34, 39):
        init[i] = np.random.uniform(0.0000001, 1 - sum(init[33:i]))

    # Set the default state to match with the input and reset
    s.set_default_state(init)
    s.reset()

    # Run to steady-state
    try:
        s.pre(2000000)
        run_sim = s.run(10, log=['intracellular_ions.ki', 'intracellular_ions.kss', 'intracellular_ions.nai', 'intracellular_ions.nass', 'membrane.V', 'intracellular_ions.cai', 'intracellular_ions.cass', 'intracellular_ions.cajsr', 'intracellular_ions.cansr'], log_times = [0])
        return run_sim
    
    except:
        return None



# Set G0 values to explore
nb_planes = 10
nb_points_per_axis = 10
nb_points = nb_points_per_axis**2

G0s = np.linspace(-13, 20, nb_planes)
samples = np.zeros((nb_planes, nb_points, 8))
states = np.zeros((nb_planes, nb_points, 9))

# For each G0 hyperplane, eavaluate the steady-state
for plane in range(nb_planes):
    print('\n Plane : ' + str(plane) + '\n')
    G0 = G0s[plane]
    
    # Set G0 in the simulation
    s.set_constant('membrane.c0', G0)
    
    # Compute the initial conditions to cover the hyperplanes
    samples[plane, :, :] = sample_hyperplane(G0, nb_points_per_axis)
    
    for sample in range(nb_points):
        # Run the simulation for each sample
        sim_res = run_sim(samples[plane, sample, :])
        
        if sim_res is not None:
            # Save the results to the array
            states[plane, sample,:] = [sim_res['intracellular_ions.nai'][-1],
                   sim_res['intracellular_ions.nass'][-1],
                   sim_res['intracellular_ions.ki'][-1],
                   sim_res['intracellular_ions.kss'][-1],
                   sim_res['intracellular_ions.cai'][-1],
                   sim_res['intracellular_ions.cass'][-1],
                   sim_res['intracellular_ions.cansr'][-1],
                   sim_res['intracellular_ions.cajsr'][-1],
                   sim_res['membrane.V'][-1]]
        else:
            print('Simulation for sample ' + str(sample) + ' didnt work.')
            
        # Log the progress
        if sample % nb_points_per_axis == 0:
            print('Sample : ' + str(sample))


# In[]
def plot_hyperplane(plane, ax):    
    # Prepare the colorscale
    viridis = cm.get_cmap('Spectral_r', 500)

    # Nai for Y
    Y_low = 0
    Y_up = 16
    Y = np.array([np.linspace(Y_low, Y_up, 100), 
                  np.linspace(Y_low, Y_up, 100)])
                
    # Cai for Z
    Z_low = 0
    Z_up = 1.5e-04
    Z = np.array([np.linspace(Z_low, Z_low, 100), 
                  np.linspace(Z_up, Z_up, 100)])
    
    G0 = G0s[plane]
    
    # Ki for X
    # Compute Ki for the extremes
    X1_min = Ki_calc(G0, Nai = Y_low, Nass = Y_low, Cai = Z_low, Cass = default[6], Cansr = default[7], Cajsr = default[8], V=-88)
    X1_max = Ki_calc(G0, Nai = Y_up, Nass = Y_up, Cai = Z_low, Cass = default[6], Cansr = default[7], Cajsr = default[8], V=-88)
    X2_min = Ki_calc(G0, Nai = Y_low, Nass = Y_low, Cai = Z_up, Cass = default[6], Cansr = default[7], Cajsr = default[8], V=-88)
    X2_max = Ki_calc(G0, Nai = Y_up, Nass = Y_up, Cai = Z_up, Cass = default[6], Cansr = default[7], Cajsr = default[8], V=-88)
    
    # Deduce the array to plot the surface
    X = np.array([np.linspace(X1_min, X1_max, 100),
                  np.linspace(X2_min, X2_max, 100)])
            
    # Take out the point outside of the plotting limits
    X[X>150] = np.NaN
    X[X<120] = np.NaN
    
    # Plot the plane for initial conditions
    ax.plot_surface(X, Y, Z, color =  col(plane), alpha = 0.1)
    
    # Set the view and limits for axes
    ax.view_init(elev=20, azim=-160)
    ax.set_xlim([120, 150])
    ax.set_ylim([Y_low, Y_up])
    ax.set_zlim([Z_low, Z_up])
    
    return None


# In[Plot the results for ORd CiPA]
# Prepare the colorscale
viridis = cm.get_cmap('Spectral_r', 500)

# Prepare the figure frame
fig = plt.figure(figsize = (13, 13))
ax = fig.add_subplot(111, projection='3d')

first_point = 0
last_point = -2

def col(index):
    return viridis(1 - (index-first_point)/(nb_planes+last_point-first_point))

for plane in range(first_point, nb_planes + last_point + 1):
    plot_hyperplane(plane, ax)
    # Plot the steady-state points in the plot range
    for sample in range(nb_points):
        ax.scatter(states[plane, sample, 2], states[plane, sample, 0], states[plane, sample, 4],
                   linewidth = 5, color = col(plane))
    
    
# Add the projections on both axes for the extremes
# Vertical
points =  np.array([[states[last_point, last_point, 2], states[last_point, last_point, 0], 0], [states[last_point, last_point, 2], states[last_point, last_point, 0], states[last_point, last_point, 4]]])
ax.plot(points[:, 0], points[:, 1], points[:, 2],
           linewidth = 3, color = col(nb_planes + last_point), LineStyle = '--')
points =  np.array([[states[first_point, 0, 2], states[first_point, 0, 0], 0], [states[first_point, 0, 2], states[first_point, 0, 0], states[first_point, 0, 4]]])
ax.plot(points[:, 0], points[:, 1], points[:, 2],
           linewidth = 3, color = col(first_point), LineStyle = '--')

# Horizontal
Ki_high = Ki_calc(G0s[last_point],
                  Nai = 0, 
                  Nass =  0, 
                  Cai =  states[last_point, last_point, 4], 
                  Cansr =  states[last_point, last_point, 5],
                  Cajsr =  states[last_point, last_point, 6], 
                  Cass =  states[last_point, last_point, 7])
Nai_high = Ki_calc(G0s[last_point],
                  Nai = 120, 
                  Nass =  120, 
                  Cai =  states[last_point, last_point, 4], 
                  Cansr =  states[last_point, last_point, 5],
                  Cajsr =  states[last_point, last_point, 6], 
                  Cass =  states[last_point, last_point, 7])

points =  np.array([[120, Nai_high, states[last_point, last_point, 4]], [Ki_high, 0, states[last_point, last_point, 4]]])

ax.plot(points[:, 0], points[:, 1], points[:, 2],
        linewidth = 3, color = col(nb_planes + last_point), LineStyle = '--')

Ki_high = Ki_calc(G0s[first_point],
                  Nai = 16, 
                  Nass =  16, 
                  Cai =  states[first_point, first_point, 4], 
                  Cansr =  states[first_point, first_point, 5],
                  Cajsr =  states[first_point, first_point, 6], 
                  Cass =  states[first_point, first_point, 7])
Nai_high = Ki_calc(G0s[first_point],
                  Nai = 150, 
                  Nass =  150, 
                  Cai =  states[first_point, first_point, 4], 
                  Cansr =  states[first_point, first_point, 5],
                  Cajsr =  states[first_point, first_point, 6], 
                  Cass =  states[first_point, first_point, 7])
points =  np.array([[150, Nai_high, states[first_point, first_point, 4]], [Ki_high, 16, states[first_point, first_point, 4]]])

ax.plot(points[:, 0], points[:, 1], points[:, 2], 
        linewidth = 3, color = col(first_point), LineStyle = '--')


# Set Ca2+ ticks
ax.set_zticks(np.linspace(0, 0.00014, 3))
ax.set_zticklabels(np.linspace(0, 0.14, 3))


ax.view_init(elev=20, azim=-160)
ax.tick_params(axis='both', labelsize= 23)
ax.set_xlabel('$[K^+]_i$ (mM)', fontsize = 35, labelpad = 25)
ax.set_ylabel('$[Na^+]_i$ (mM)', fontsize = 35, labelpad = 25)
ax.set_zlabel(r'$[Ca^{2+}]_i$ ($\mu$M)', fontsize = 35, labelpad = 25)

ax.set_xlim([120, 150])
ax.set_ylim([0, 16])
ax.set_zlim([0, 1.5e-04])


