# -*- coding: utf-8 -*-
"""
Created on Tue Jan 21 14:25:17 2020

@author: yann-stanislas.barral@roche.com
"""


import sabs_pkpd
import matplotlib.pyplot as plt
import myokit
import numpy as np


# In[Load model]
filename = 'C:/Users/barraly/Documents/PhD/Sensitivity to initial conditions/OHara CiPA.mmt'
s = sabs_pkpd.load_model.load_simulation_from_mmt(filename)
s.set_tolerance(abs_tol = 1e-08, rel_tol = 1e-06)
default_state = s.state()


# In[Select initial conditions]
# Set 1 Hz pacing protocol
p = myokit.Protocol()
e = myokit.ProtocolEvent(level=1, start=50, duration=1, period=1000, multiplier=0)
p.add(e)
s.set_protocol(p)

# Reset the model
s.reset()
s.set_state(default_state)

# Run the model to 1 Hz steady-state
s.pre(1200000)
one_hz_steady = s.run(100000, log_interval = 10)
one_hz = s.run(1000, log_interval = 1)
one_hz_norm_Ki = (one_hz['intracellular_ions.ki'] - np.min(one_hz['intracellular_ions.ki'])) / (np.max(one_hz['intracellular_ions.ki']) - np.min(one_hz['intracellular_ions.ki']))
new_default_state = s.state()

# Add perturbation: heart rate to 2 Hz
p = myokit.Protocol()
e = myokit.ProtocolEvent(level=1, start=50, duration=1, period=500, multiplier=0)
p.add(e)
s.set_protocol(p)

# Observe adaptation to HR
s.reset()
s.set_state(new_default_state)
adapt_HR = s.run(1000000, log_interval = 10)

# Retrieve the late convergence with finer resolution for the phase portrait
s.reset()
s.set_state(new_default_state)
s.pre(500000)
adapt_HR_phase = s.run(600000, log_interval = 1)
two_hz = s.run(1000, log_interval = 1)
two_hz_norm_Ki = (two_hz['intracellular_ions.ki'] - np.min(two_hz['intracellular_ions.ki'])) / (np.max(two_hz['intracellular_ions.ki']) - np.min(two_hz['intracellular_ions.ki']))

#
## Add perturbation: hypokalemia
#s.set_constant('extracellular.ko', 3.5)
#a = s.run(2000000, log_interval = 1000)
#hypokalemia_convergence = a['intracellular_ions.ki']


# In[Plot the results]
f, (a0, a1, a2) = plt.subplots(3, 2, gridspec_kw={'height_ratios': [1, 1, 2.5]})
f.set_figheight(15)
f.set_figwidth(15)

# Plot the voltage profile at steady-state
a0[0].plot(np.linspace(0, 0.999, 1000), one_hz['membrane.V'], label = '1 Hz', LineWidth = 4, color = 'k')
a0[0].plot(np.linspace(0, 0.999, 1000), two_hz['membrane.V'], label = '2 Hz', LineWidth = 4, color = 'b')
a0[0].legend(fontsize = 22, loc = 'upper right')
a0[0].set_ylabel('Voltage (mV)', Fontsize = 25)
a0[0].set_yticks([-80, -40, 0, 40])
a0[0].tick_params(axis='both', which='major', labelsize=20)

# Plot the K+ profile at steady-state (relative)
a1[0].plot(np.linspace(0, 0.999, 1000), one_hz_norm_Ki, label = '1 Hz', LineWidth = 4, color = 'k')
a1[0].plot(np.linspace(0, 0.999, 1000), two_hz_norm_Ki, label = '2 Hz', LineWidth = 4, color = 'b')
#a1[0].hlines(y = one_hz_norm_Ki[0], xmin = 0, xmax = 0.999, LineWidth = 2, alpha = 0.2, color = 'k', LineStyle = '--')
#a1[0].hlines(y = two_hz_norm_Ki[0], xmin = 0, xmax = 0.999, LineWidth = 2, alpha = 0.2, color = 'b', LineStyle = '--')
a1[0].scatter(0.05, one_hz_norm_Ki[50], color = 'k', s = 100, marker='o', LineWidth = 2.5)
a1[0].scatter([0.05, 0.55], [two_hz_norm_Ki[50], two_hz_norm_Ki[550]], color = 'b', s = 100, marker='o', LineWidth = 2.5)
a1[0].legend(fontsize = 22)
a1[0].set_ylabel('Relative $[K^+]_i$', Fontsize = 23)
a1[0].tick_params(axis='both', which='major', labelsize=20)
a1[0].set_xlabel('Time (s)', Fontsize = 25)
a1[0].set_yticks([0, 1])
a1[0].set_yticklabels(['min', 'max'])

# Plot the adapation of V to HR
a0[1].plot(np.linspace(96, 99.990, 400), one_hz_steady['membrane.V'][:400], color = 'k', LineWidth = 3)
a0[1].plot(np.linspace(100, 105.990, 600),adapt_HR['membrane.V'][:600], color = 'r', LineWidth = 3)
a0[1].scatter(np.linspace(96, 99, 4), one_hz_steady['membrane.V'][:400:100], color = 'k', s = 100, marker='o', LineWidth = 2.5)
a0[1].scatter(np.linspace(100, 105.5, 12), adapt_HR['membrane.V'][:600:50], color = 'r', s = 100, marker='o', LineWidth = 2.5)
a0[1].set_ylabel('Voltage (mV)', Fontsize = 23)
#a0[1].set_xlabel('Time (s)', Fontsize = 23)
a0[1].tick_params(axis='both', which='major', labelsize=20)

# Plot the adapation of K+ to HR
a1[1].plot(np.linspace(96, 99.990, 400), one_hz_steady['intracellular_ions.ki'][:400], color = 'k', LineWidth = 3)
a1[1].plot(np.linspace(100, 105.990, 600), adapt_HR['intracellular_ions.ki'][:600], color = 'r', LineWidth = 3)
a1[1].scatter(np.linspace(96, 99, 4), one_hz_steady['intracellular_ions.ki'][:400:100], color = 'k', s = 100, marker='o', LineWidth = 2.5)
a1[1].scatter(np.linspace(100, 105.5, 12), adapt_HR['intracellular_ions.ki'][:600:50], color = 'r', s = 100, marker='o', LineWidth = 2.5)
a1[1].set_ylabel('$[K^+]_i$ (mM)', Fontsize = 23)
a1[1].set_xlabel('Time (s)', Fontsize = 23)
a1[1].ticklabel_format(useOffset=False)
a1[1].tick_params(axis='both', which='major', labelsize=20)


# Plot the phase portrait of V against Ki
a2[0].plot(adapt_HR_phase['intracellular_ions.ki'][:], adapt_HR_phase['membrane.V'][:], color = 'r', LineWidth = 3, label = 'Transient')
#a2[0].scatter(adapt_HR['intracellular_ions.ki'][50000::50], adapt_HR['membrane.V'][50000::50], color = 'r', s = 100, marker='o', LineWidth = 2.5)
#a2[0].plot(one_hz_steady['intracellular_ions.ki'][:], one_hz_steady['membrane.V'][:], color = 'k', LineWidth = 3, label = '1 Hz steady-state')
#a2[0].scatter(one_hz_steady['intracellular_ions.ki'][:100:100], one_hz_steady['membrane.V'][:100:100], color = 'k', s = 100, marker='o', LineWidth = 2.5)
a2[0].plot(two_hz['intracellular_ions.ki'], two_hz['membrane.V'], color = 'b', LineWidth = 3, label = '2 Hz steady-state')
a2[0].scatter(two_hz['intracellular_ions.ki'][50], two_hz['membrane.V'][50], color = 'b', s = 100, marker='o', LineWidth = 2.5)
a2[0].set_ylabel('Voltage (mV)', Fontsize = 23)
a2[0].set_xlabel('$[K^+]_i$ (mM)', Fontsize = 23)
a2[0].tick_params(axis='both', which='major', labelsize=20)
a2[0].set_xticks([143.276, 143.280, 143.284, 143.288])
a2[0].set_yticks([-80, -40, 0, 40])
a2[0].ticklabel_format(useOffset=False)
a2[0].legend(fontsize = 20, loc = 'upper left')


# Plot the full of convergence of K+ to new steady-state
a2[1].plot(np.linspace(0, 99, 100), one_hz_steady['intracellular_ions.ki'][::100], color = 'k', LineWidth = 3)
a2[1].plot(np.linspace(100, 799.5, 1200), adapt_HR['intracellular_ions.ki'][:60000:50], color = 'r', LineWidth = 3)
a2[1].plot(np.linspace(800, 899.5, 200), adapt_HR['intracellular_ions.ki'][90000:100000:50], color = 'b', LineWidth = 3)
a2[1].set_ylabel('Diastolic $[K^+]_i$ (mM)', Fontsize = 23)
a2[1].set_xlabel('Time (s)', Fontsize = 23)
a2[1].ticklabel_format(useOffset=False)
a2[1].tick_params(axis='both', which='major', labelsize=20)
a2[1].vlines(ymin = 143.3, ymax = 144.5, x = 100, LineStyle = '--', color = 'k', LineWidth = 3)
#a2[1].vlines(ymin = 143.3, ymax = 144.5, x = 600, LineStyle = '--', color = 'k', LineWidth = 3)
a2[1].text(x = -25, y = 143.8, s = '1 Hz', Fontsize = 20)
a2[1].text(x = 300, y = 143.8, s = '2 Hz', Fontsize = 20)

plt.tight_layout()

