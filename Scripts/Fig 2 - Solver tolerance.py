# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 14:25:17 2020

@author: yann-stanislas.barral@roche.com
"""


import sabs_pkpd
import matplotlib.pyplot as plt
import os
import numpy as np

# Select the folder in which this repo is downloaded in the line below
os.chdir('../')


# In[Load model]
filename = './Models/OHara CiPA.mmt'
s = sabs_pkpd.load_model.load_simulation_from_mmt(filename)

runtime = 2000000
log_interval = 1000

# In[Select initial conditions]
# Names of the state variables in the MMT model
# Note that V is added in the list here to be called as a read out

# Set the solver tolerance to high
s.set_tolerance(abs_tol = 1e-06, rel_tol = 1e-04)
# Reset the model to the initial conditions retrieved when loading the model
s.set_state(sabs_pkpd.constants.default_state)
# run the model to limit cycle
a = s.run(runtime, log_interval = log_interval)


# Set the solver tolerance to low
s.set_tolerance(abs_tol = 1e-08, rel_tol = 1e-06)
# Reset the model to the initial conditions retrieved when loading the model
s.set_state(sabs_pkpd.constants.default_state)
# run the model to limit cycle
a1 = s.run(runtime, log_interval = log_interval)


# In[Plot the results]
state_name = 'membrane.V'

offset = 250
y = a[state_name][offset:]
y1 = a1[state_name][offset:]
x = np.arange(offset, offset + len(y))

plt.figure(figsize = (7.5, 7.5))
plt.plot(x, y, label = 'abs_tol = 1e-06, rel_tol = 1e-04', linewidth = 5)
plt.plot(x, y1, label = 'abs_tol = 1e-08, rel_tol = 1e-06', linewidth = 5)
plt.xlabel('# of paces', fontsize = 15)
plt.ylabel('Resting membrane potential (mV)', fontsize = 15)
plt.xticks(fontsize = 15)
plt.legend(fontsize = 15)

# Save
plt.tight_layout()
plt.savefig('./Figures/Numerical error.png')

