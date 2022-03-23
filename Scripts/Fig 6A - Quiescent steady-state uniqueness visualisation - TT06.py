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
s = sabs_pkpd.load_model.load_simulation_from_mmt('./Models/tentusscher_2006 - analytical voltage.mmt')
default = s.state()

# Set solver tolerance to fine
s.set_tolerance(abs_tol=1e-08, rel_tol = 1e-08)

# Set no pacing for the full model
s.set_constant('stimulus.amplitude', 0)


# In[For TT06]
def G0_calc(Ki, Nai, Cai, Casr, Cass, V=-85.3, extraK=5.4, extraNa=140, extraCa=2):    
    tot_cai = Cai * (1 + 0.2 / (Cai + 0.001))
    tot_cass = Cass * (1 + 0.4 / (Cass + 0.00025))
    tot_casr = Casr * (1 + 10 / (Casr + 0.3))
    Vmyo = 16.404
    Vsr = 1.094
    Vss = 5.468e-2
    tot_ca = tot_cai + tot_cass * Vss/Vmyo + tot_casr * Vsr/Vmyo
    G0 = V / (96.4853415 * Vmyo) * 0.185 - (Ki + Nai + 2 * tot_ca - extraK - extraNa - 2 * extraCa)
    return G0

# Define the function to sample a G0 plane
def sample_hyperplane(G0, nb_points_per_axis = 15):
    nb_points = nb_points_per_axis**2
    points = np.zeros((nb_points, 5))
    
    for i in range(nb_points_per_axis):
        # Cai
        points[nb_points_per_axis * i:nb_points_per_axis * (i+1), 0] = np.linspace(0.5, 1.5, nb_points_per_axis) * default[0]
        #Casr
        points[nb_points_per_axis * i:nb_points_per_axis * (i+1), 1] = np.random.uniform(0.5, 1.5, size = nb_points_per_axis) * default[1]
        #Cass
        points[nb_points_per_axis * i:nb_points_per_axis * (i+1), 2] = np.random.uniform(0.5, 1.5, size = nb_points_per_axis) * default[2]
        # Nais to be sampled from the physiological range of concentrations
        points[nb_points_per_axis * i:nb_points_per_axis * (i+1), 3] = np.linspace(4, 16, nb_points_per_axis)
    
    # Compute Ki for each point now to stick to the same G0 value
    for sample in range(nb_points):
        points[sample, 4] = Ki_calc(G0,
              Cai = points[sample, 0],                     
              Casr = points[sample, 1],
              Cass = points[sample, 2],
              Nai = points[sample, 3],  
              V=-85.3)
    return points

# Define the function to run the model to steady-state from 
# a given starting point
def run_sim(concs):
    # Set the default state to match with the input and reset
    init = default.copy()
    init[0] = concs[0]
    init[1] = concs[1]
    init[2] = concs[2]
    init[4] = concs[3]
    init[5] = concs[4]
    
    # Set non-concentrations state variables to random
    init[3] = np.random.uniform(0.0000001, 0.999999)
    init[6:] = np.random.uniform(0.0000001, 0.999999, size = 12)
    
    s.set_default_state(init)
    s.reset()

    # Run to steady-state wihtout pacing
    run_sim = s.run(4000001, log=['potassium.K_i', 'sodium.Na_i', 'membrane.V', 'calcium.Ca_i', 'calcium.Ca_sr', 'calcium.Ca_ss'], log_times = [0, 4000000])
    
    return run_sim


# In[Execute the functions for TT06]
# Define how to compute Ki for TT06    
def Ki_calc(G0, Nai, Cai, Casr, Cass, V = -85.3, extraK = 5.4, extraNa = 140,
            extraCa = 2):
    Cai_tot = Cai + Cai*0.2 / (Cai + 0.001)
    Casr_tot = Casr + Casr * 10 / (Casr + 0.3)
    Cass_tot = Cass + Cass * 0.4 / (Cass + 0.00025)
    Ki = extraK + extraNa + 2 * extraCa - G0 + V / 96.4853415 * 0.185 / 16.404 - (Nai + 2*Cai_tot  + 2*0.001094/0.016404*Casr_tot + 2*5.468e-5/0.016404*Cass_tot)
    return Ki

# Set G0 values to explore
nb_planes = 10
nb_points_per_axis = 10
nb_points = nb_points_per_axis**2

G0s = np.linspace(-20.8, 13.2, nb_planes)
samples = np.zeros((nb_planes, nb_points, 5))
states = np.zeros((nb_planes, nb_points, 6))

# For each G0 hyperplane, eavaluate the quiescent steady-state
for plane in range(nb_planes):
    print('\n Plane :' + str(plane) + '\n')
    G0 = G0s[plane]
    
    # Set G0 in the simulation
    s.set_constant('membrane.c0', G0)
    
    # Compute the initial conditions to cover the hyperplanes
    samples[plane, :, :] = sample_hyperplane(G0, nb_points_per_axis)
    
    for sample in range(nb_points):
        # Run the simulation for each sample
        sim_res = run_sim(samples[plane, sample, :])
        
        # Save the results to the array
        states[plane, sample,:] = [sim_res['calcium.Ca_i'][-1],
               sim_res['calcium.Ca_sr'][-1],
               sim_res['calcium.Ca_ss'][-1],
               sim_res['sodium.Na_i'][-1],
               sim_res['potassium.K_i'][-1],
               sim_res['membrane.V'][-1]]
         
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
    
    # Ki for X
    # Compute Ki for the extremes
    X1_min = Ki_calc(G0s[plane], Cai = Z_low, Nai = Y_low, Casr = default[1], Cass = default[2])
    X1_max = Ki_calc(G0s[plane], Cai = Z_low, Nai = Y_up, Casr = default[1], Cass = default[2])
    X2_min = Ki_calc(G0s[plane], Cai = Z_up, Nai = Y_low, Casr = default[1], Cass = default[2])
    X2_max = Ki_calc(G0s[plane], Cai = Z_up, Nai = Y_up, Casr = default[1], Cass = default[2])
    
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


# In[Plot the results for TT06
# Prepare the colorscale
    
viridis = cm.get_cmap('Spectral_r', 500)

# Prepare the figure frame
fig = plt.figure(figsize = (13, 13))
ax = fig.add_subplot(111, projection='3d')

first_point = 4
last_point = -1

def col(index):
    return viridis(1 - (index-first_point)/(nb_planes+last_point - first_point))

for plane in range(first_point, nb_planes + last_point + 1):
    plot_hyperplane(plane, ax)
    # Plot the steady-state points in the plot range
    for sample in range(nb_points):
        ax.scatter(states[plane, sample, 4], states[plane, sample, 3], states[plane, sample, 0],
                   linewidth = 5, color = col(plane))
    
    
# Add the projections on both axes for the extremes
# Vertical
Ki_proj = Ki_calc(G0s[last_point], states[last_point, last_point, 3], 0, states[last_point, last_point, 1], states[last_point, last_point, 2])
points =  np.array([[Ki_proj, states[last_point, last_point, 3], 0], [Ki_proj, states[last_point, last_point, 3], states[last_point, last_point, 0]]])
ax.plot(points[:, 0], points[:, 1], points[:, 2],
        linewidth = 3, color = col(nb_planes + last_point), LineStyle = '--')
Ki_proj = Ki_calc(G0s[first_point], states[first_point, first_point, 3], 0, states[first_point, first_point, 1], states[first_point, first_point, 2])
points =  np.array([[Ki_proj, states[first_point, 0, 3], 0], [Ki_proj, states[first_point, 0, 3], states[first_point, 0, 0]]])
ax.plot(points[:, 0], points[:, 1], points[:, 2],
        linewidth = 3, color = col(first_point), LineStyle = '--')

# Horizontal
Ki_high = Ki_calc(G0s[last_point], 0, states[last_point, last_point, 0], states[last_point, last_point, 1], states[last_point, last_point, 2])
Nai_high = Ki_calc(G0s[last_point], 120, states[last_point, last_point, 0], states[last_point, last_point, 1], states[last_point, last_point, 2])
points =  np.array([[120, Nai_high, states[last_point, last_point, 0]], [Ki_high, 0, states[last_point, last_point, 0]]])
ax.plot(points[:, 0], points[:, 1], points[:, 2],
        linewidth = 3, color = col(nb_planes + last_point), LineStyle = '--')

Ki_high = Ki_calc(G0s[first_point], 16, states[first_point, first_point, 0], states[first_point, first_point, 1], states[first_point, first_point, 2])
Nai_high = Ki_calc(G0s[first_point], 150, states[first_point, first_point, 0], states[first_point, first_point, 1], states[first_point, first_point, 2])
points =  np.array([[150, Nai_high, states[first_point, first_point, 0]], [Ki_high, 16, states[first_point, first_point, 0]]])
ax.plot(points[:, 0], points[:, 1], points[:, 2], 
        linewidth = 3, color = col(first_point), LineStyle = '--')

# Tick the Ca2+ axis in uM
ax.set_zticks(np.linspace(0, 0.00014, 3))
ax.set_zticklabels(np.linspace(0, 0.14, 3))

ax.tick_params(axis='both', labelsize= 23)
ax.set_xlabel('$[K^+]_i$ (mM)', fontsize = 35, labelpad = 25)
ax.set_ylabel('$[Na^+]_i$ (mM)', fontsize = 35, labelpad = 25)
ax.set_zlabel(r'$[Ca^{2+}]_i$ ($\mu$M)', fontsize = 35, labelpad = 25)
ax.view_init(elev=20, azim=-160)

#ax.ticklabel_format(style='sci', axis='z', scilimits=(0,0), useOffset=False)
ax.set_xlim([120, 150])
ax.set_ylim([0, 16])
ax.set_zlim([0, 1.5e-04])




