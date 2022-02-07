# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 14:25:17 2020

@author: barraly
"""


import sabs_pkpd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as clrs
import matplotlib.cm as cm
import matplotlib.axes
import myokit
import os

# Select the folder in which this repo is downloaded in the line below
os.chdir('The/location/of/the/root/folder/of/this/repo')


# In[Load model]
filename = './Models/Ohara CiPA - analytical voltage.mmt'
s = sabs_pkpd.load_model.load_simulation_from_mmt(filename)
default_state = s.state()
 

# In[Verify that the model is run to steady-state]
def difference(state1, state2):
    a = np.array(state1)
    b = np.array(state2)
    return np.sum(np.square(a-b))

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

def Nai_calc(G0, Ki = 144.65559, Kss = 144.65556, Cai = 8.6e-5, Cansr = 1.61957, Cajsr = 1.571234014, Cass = 8.49e-5, V=-88, extraK = 5.4, extraNa = 140, extraCa = 1.8):
    tot_cai = Cai * (1 + 0.05 / (Cai + 0.00238) + 0.07/(Cai + 0.0005))
    tot_cass = Cass * (1 + 0.047 / (Cass + 0.00087) + 1.124/(Cass + 0.0087))
    tot_cajsr = Cajsr * (1 + 10 / (Cajsr + 0.8))
    return (V / (96485 * 2.583592e-05) * 0.0001533576 + extraK + extraNa + 2 * extraCa - G0 - Ki - Kss * 0.029411764705882353 - 2*(tot_cai + tot_cass * 0.029411764705882353 + Cansr * 0.08117647059 + tot_cajsr * 0.007059)) / 1.029411764705882353 

def compute(Ki, Nai):
    # Reinitialise the myokit.Simulation
    s.reset()
    
    # Set the initial conditions
    initial_state = default_state.copy()
    Gamma_0 = G0_calc(Ki, Ki, Nai, Nai)
    initial_state[1] = Nai
    initial_state[2] = Nai
    initial_state[3] = Ki
    initial_state[4] = Ki
    s.set_state(initial_state)
    
    # Set the value of Gamma_0 in the myokit.Simulation
    s.set_constant('membrane.c0', Gamma_0)
    
    # Record the action potential at the limit cycle
    s.pre(1500000)
    out = s.run(1000, log_interval = 1)
    states_end_1 = s.state()
    out_2 = s.run(1000, log_interval = 1)
    states_end_2 = s.state()
    
    # In case of alternans, select the shorter AP
    APD90_1 = sabs_pkpd.cardiac.compute_APD(np.array(out['membrane.V']), time_points=np.linspace(0, 999, 1000), upstroke_time=50)
    APD90_2 = sabs_pkpd.cardiac.compute_APD(np.array(out_2['membrane.V']), time_points=np.linspace(0, 999, 1000), upstroke_time=50)
    
    if APD90_1 < APD90_2:
        AP = np.array(out['membrane.V'])
        Cai = out['intracellular_ions.cai']
        states_end = states_end_1
        
    else:
        AP = np.array(out_2['membrane.V'])
        Cai = out_2['intracellular_ions.cai']
        states_end = states_end_2
        
    return states_end, AP, Cai


# In[Initialise the arrays]
# Number of C0 lines
number_of_lines = 15

# Number of points taken per C0 line
number_init = 10

# Literature boundaries for intracellular potassium and sodium
Ki_min = 120
Ki_max = 152
Nai_min = 4
Nai_max = 16

# Compute corresponding range of C0 that verify V_initial = -88 mV
G0s = np.linspace(int(-22), int(20), number_of_lines)

# Compute the arrays of (Ki, Nai) to sample the square boundaries for initial conditions
Kis = np.zeros((number_of_lines, number_init))
Nais = np.zeros((number_of_lines, number_init))

for i in range(number_of_lines):
    # Compute the limits of each line that remains within the physiological area
    low_lim_nai = Ki_calc(G0s[i], Nai_max, Nai_max, default_state[5], default_state[7], default_state[8], default_state[6])
    low_Ki = max(Ki_min, low_lim_nai)
    up_Nai = Nai_calc(G0s[i], low_Ki, low_Ki, default_state[5], default_state[7], default_state[8], default_state[6])
    
    up_lim_nai = Ki_calc(G0s[i], Nai_min, Nai_min, default_state[5], default_state[7], default_state[8], default_state[6])
    up_Ki = min(Ki_max, up_lim_nai )
    low_Nai = Nai_calc(G0s[i], up_Ki, up_Ki, default_state[5], default_state[7], default_state[8], default_state[6])
    
    Kis[i,:] = np.linspace(low_Ki, up_Ki, number_init)
    Nais[i,:] = np.linspace(up_Nai, low_Nai, number_init)
    
# Visualise the lines of sampling
plt.figure(figsize = (7, 5))
plt.xlabel('Nai (mM)', fontsize = 22)
plt.ylabel('Ki (mM)', fontsize = 22)
for i in range(number_of_lines):
    plt.scatter(Nais[i, :], Kis[i, :], color = 'k')


# In[Run the model]
def initialise_arrays():
    # Arrays to record simulations outputs
    inits = np.zeros((number_of_lines, number_init, len(default_state)))
    scores = np.zeros((number_of_lines, number_init))
    states = np.zeros((number_of_lines, number_init, len(default_state)))
    APs = np.zeros((number_of_lines, number_init, 1000))
    APDs = []
    Cais = np.zeros((number_of_lines, number_init, 1000))
    CaTDs = []
    
    return inits, scores, states, APs, APDs, Cais, CaTDs


def Nai_Ki_map(hERG_rescale):
    inits, scores, states, APs, APDs, Cais, CaTDs = initialise_arrays()
    
    # Set the drug block effect in the model
    s.set_constant('drug.ikr_rescale', hERG_rescale)
    
    for i in range(number_of_lines):
        #Run the model to steady-state, starting from points sampled on each C0 line
        print('Computing line ' + str(i))
        
        for j in range(number_init):
            # Log the progress
            print('j = ' + str(j))
            
            # Run the model to steady-state
            states_end, AP, Cai = compute(Kis[i, j], Nais[i, j])
            
            # Save the outputs to the corresponding arrays
            states[i, j, :] = states_end
            APs[i, j, :] = AP
            Cais[i, j, :] = Cai
            
            # Compute the metrics
            APD = sabs_pkpd.cardiac.compute_APD(APs[i, j, :], time_points = np.linspace(0, 999, 1000), upstroke_time = 50, print_warnings=False)
            APDs.append(APD)
            CaTD = sabs_pkpd.cardiac.compute_calcium_transient_duration(Cais[i, j, :], time_points = np.linspace(0, 999, 1000), upstroke_time = 50, print_warnings=False)
            CaTDs.append(CaTD)

    return inits, scores, states, APs, APDs, Cais, CaTDs


# In[]
inits_all = []
scores_all = []
states_all = []
APs_all = []
APDs_all = []
Cais_all = []
CaTDs_all = []

# Try various hERG rescales
for hERG_rescale in [1, 0.8, 0.6, 0.4, 0.2, 0.1, 0.05]:
    inits, scores, states, APs, APDs, Cais, CaTDs = Nai_Ki_map(hERG_rescale)
    plot_2D_Nai_Ki_map(hERG_rescale, states)
    
    # Save to a bigger structure
    inits_all.append(inits)
    scores_all.append(scores)
    states_all.append(states)
    APs_all.append(APs)
    APDs_all.append(APDs)
    Cais_all.append(Cais)
    CaTDs_all.append(CaTDs)
    

# In[Plot Nai Ki map for initial lines and steady-states reached]
def place_caption_label(ax, label, loc='upper left', fontsize=35):
    from matplotlib.offsetbox import AnchoredText
    at = AnchoredText(label, loc=loc, prop=dict(size=fontsize), frameon=True, borderpad = 0)
    ax.add_artist(at)    
    return None

def plot_2D_Nai_Ki_map(hERG_rescale, states, label_each_line = False, ax=None):
    if ax is None:
        fig, ax = plt.figure(figsize=(15,15))
        ax.title('Nai / Ki map at initial state, ' + str(int((1-hERG_rescale)*100)) + '% of hERG block', fontsize = 30)
        
    ax.set_xlabel('$[Na^+]_i$ (mM)', fontsize = 60)
    ax.set_ylabel('$[K^+]_i$ (mM)', fontsize = 60)
    ax.set_xlim([4, 16])
    ax.set_ylim([120, 152])
    
    ax.tick_params(labelsize = 40)
    
    cmap = cm.get_cmap('gnuplot', number_of_lines)
    newcolors = cmap(np.linspace(0, 1, number_of_lines))
    cmap = clrs.ListedColormap(newcolors)
    
    ax.hlines(Ki_min, xmin = (Nai_min-4)/16, xmax = (Nai_max-4)/16, color = 'k', linestyle = '--')
    ax.hlines(Ki_max, xmin = (Nai_min-4)/16, xmax = (Nai_max-4)/16, color = 'k', linestyle = '--')
    ax.vlines(Nai_min, ymin = (Ki_min-110)/50, ymax = (Ki_max-110)/50, color = 'k', linestyle = '--')
    ax.vlines(Nai_max, ymin = (Ki_min-110)/50, ymax = (Ki_max-110)/50, color = 'k', linestyle = '--')
    
    ax.plot([Nais[0,0], Nais[0, -1]], [Kis[0, 0], Kis[0, -1]], color = cmap(0), linewidth = 3, label = 'Initial conditions')
    ax.scatter(states[0, :, 2], states[0, :, 4], edgecolors = cmap(0), linewidth = 2, s = 150, facecolors='none', label = 'Steady-state')
    
    for i in range(number_of_lines):
        ax.plot([Nais[i,0], Nais[i, -1]], [Kis[i, 0], Kis[i, -1]], color = cmap(number_of_lines- 1 - i), linewidth = 3)
        ax.scatter(states[i, :, 2], states[i, :, 4], edgecolors = cmap(number_of_lines - 1 - i), linewidth = 2, s = 150, facecolors='none')
        
    # Label each G0 line
    if label_each_line:
        rota = -18
        size = 21
        ax.text(x = 4, y = 119.4, s = '$\Gamma_0 = $' + str(int(G0s[-1])) + ' mM', fontsize = size, rotation = rota)
        ax.text(x = 4, y = 122.4, s = '$\Gamma_0 = $' + str(int(G0s[-2])) + ' mM', fontsize = size, rotation = rota)
        ax.text(x = 9.6, y = 120, s = '$\Gamma_0 = $' + str(int(G0s[-3])) + ' mM', fontsize = size, rotation = rota)
        ax.text(x = 12.1, y = 120, s = '$\Gamma_0 = $' + str(int(G0s[-4])) + ' mM', fontsize = size, rotation = rota)
        for i in range(4, 15):
            ax.text(x = 14, y = 110 + 2.88*i, s = '$\Gamma_0 = $' + str(int(G0s[14-i])) + ' mM', fontsize = size, rotation = rota)
        
    
    # Save if needed
    if ax is None:
        ax.legend(fontsize = 25)
        folder = 'C:/Users/barraly/Documents/PhD/Sensitivity to initial conditions/Nai Ki map evolution/ORd CiPA/'
        save_filename = folder + 'Nai Ki map at steady-state, hERG rescale = ' + str(hERG_rescale) + '.png'
        plt.savefig(save_filename)

    return None

fig, axes = plt.subplots(1, 3, figsize = (30, 10))
plot_2D_Nai_Ki_map(1, states_all[0], label_each_line = True, ax=axes[0])
plot_2D_Nai_Ki_map(0.1, states_all[-2], ax=axes[1])
plot_2D_Nai_Ki_map(0.05, states_all[-1], ax=axes[2])

# Labels
axes[2].legend(fontsize = 28, loc = 'lower right')
place_caption_label(axes[0], 'A', 'upper left', 40)
place_caption_label(axes[1], 'B', 'upper left', 40)
place_caption_label(axes[2], 'C', 'upper left', 40)

# Save
plt.tight_layout()
plt.savefig('./Figures/bifurcation illustration.png', dpi = 300)

    
    