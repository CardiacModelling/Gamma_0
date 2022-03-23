# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 09:13:10 2021

@author: barraly
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import os 

# Select the folder in which this repo is downloaded in the line below
os.chdir('The/location/of/the/root/folder/of/this/repo')
    

# In[Load the data]
data = pd.read_csv('./Scripts/conc_summary.txt', sep = ' ')


# In[Plot the data]
size_labels = 28
size_ticks = 25

plt.figure(figsize = (12, 8))

#plt.subplot(1, 2, 1)
#plt.title('Initial conditions', Fontsize = 30)
plt.xlabel('$[Na^+]_i$ (mM)', fontsize = size_labels)
plt.ylabel('$[K^+]_i$ (mM)', fontsize = size_labels)
plt.xticks(fontsize = size_labels)
plt.yticks(fontsize = size_labels)
for row in range(data['Model'].count()):
    if data['Species'][row] == 'Human':
        color = '#7fc97f'
    if data['Species'][row] == 'Dog':
        color = '#beaed4'
    if data['Species'][row] == 'Rabbit':
        color = '#fdc086'
    if data['Species'][row] == 'GuineaPig':
        color = '#ffff99'
    if data['Species'][row] == 'Mammalian':
        color = '#386cb0'
    if data['Species'][row] == 'Rat':
        color = '#f0027f'
    if data['Species'][row] == 'Mouse':
        color = '#f0027f'


    if data['Tissue'][row] == 'Ventricle':
        shape = 'o'
    if data['Tissue'][row] == 'Atrium':
        shape = '+'
    if data['Tissue'][row] == 'Purkinje':
        shape = 'v'
    if data['Tissue'][row] == 'SAN':
        shape = 'x'  
     
    plt.scatter(data['Sodium_i_IC'][row], data['Potassium_i_IC'][row], color = color, marker = shape, s = 200)

# Add the delimitation of the extreme concentrations
plt.hlines(y = 120, xmin = 4, xmax = 16, linestyle = '--', color = 'k', linewidth = 4)
plt.hlines(y = 152.4, xmin = 4, xmax = 16, linestyle = '--', color = 'k', linewidth = 4)
plt.vlines(x = 4, ymin = 120, ymax = 152.4, linestyle = '--', color = 'k', linewidth = 4)
plt.vlines(x = 16, ymin = 120, ymax = 152.4, linestyle = '--', color = 'k', linewidth = 4)

# Label Tomek and Grandi models
plt.text(x= 4.1, y = 113, s = 'Grandi et al. 2010', fontsize = 21)
plt.text(x= 13.7, y = 156, s = 'Tomek et al. 2020', fontsize = 21)
plt.ylim(95, 160)

# Add the legend
plt.scatter([], [], color = 'k', marker = 'o', s = 200, label = 'Ventricle')
plt.scatter([], [], color = 'k', marker = '+', s = 200, label = 'Atrium')
plt.scatter([], [], color = 'k', marker = 'v', s = 200, label = 'Purkinje')
plt.scatter([], [], color = 'k', marker = 'x', s = 200, label = 'Sino-atrial node')
plt.legend(fontsize = 22, loc='lower right')

# Save
plt.tight_layout()
plt.savefig('./Figures/Range of concentrations.png', dpi = 300)

