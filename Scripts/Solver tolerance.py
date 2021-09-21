# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 14:25:17 2020

@author: yann-stanislas.barral@roche.com
"""


import sabs_pkpd
import matplotlib.pyplot as plt


# In[Load model]
filename = 'initial_conditions_study/Models/OHara CiPA.mmt'
s = sabs_pkpd.load_model.load_simulation_from_mmt(filename)


# In[Select initial conditions]
# Names of the state variables in the MMT model
# Note that V is added in the list here to be called as a read out

# Set the solver tolerance to high
s.set_tolerance(abs_tol = 1e-06, rel_tol = 1e-04)
# Reset the model to the initial conditions retrieved when loading the model
s.set_state(sabs_pkpd.constants.default_state)
# run the model to limit cycle
a = s.run(2000000, log_interval = 1000)


# Set the solver tolerance to low
s.set_tolerance(abs_tol = 1e-08, rel_tol = 1e-06)
# Reset the model to the initial conditions retrieved when loading the model
s.set_state(sabs_pkpd.constants.default_state)
# run the model to limit cycle
a1 = s.run(2000000, log_interval = 1000)


# In[Plot the results]
state_name = 'membrane.V'

plt.figure(figsize = (7.5, 7.5))
plt.plot(a[state_name][250:], label = 'abs_tol = 1e-06, rel_tol = 1e-04', LineWidth = 5)
plt.plot(a1[state_name][250:], label = 'abs_tol = 1e-08, rel_tol = 1e-06', LineWidth = 5)
plt.xlabel('# of paces', Fontsize = 15)
plt.ylabel('Resting membrane potential (mV)', Fontsize = 15)
plt.xticks(Fontsize = 15)
plt.legend(fontsize = 15)

