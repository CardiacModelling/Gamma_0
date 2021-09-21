# -*- coding: utf-8 -*-
"""
Created on Thu Jul  1 09:13:10 2021

@author: barraly
"""

import numpy as np
import matplotlib.pyplot as plt
import pandas as pd


# In[Load the data]
folder = 'C:/Users/barraly/Documents/PhD/Sensitivity to initial conditions/Paper elaboration/Scripts/'
data = pd.read_csv(folder + 'conc_summary.txt', sep = ' ')


# In[Plot the data]
size_labels = 28
size_ticks = 25

plt.figure(figsize = (12, 8))

#plt.subplot(1, 2, 1)
#plt.title('Initial conditions', Fontsize = 30)
plt.xlabel('$[Na^+]_i$ (mM)', Fontsize = size_labels)
plt.ylabel('$[K^+]_i$ (mM)', Fontsize = size_labels)
plt.xticks(Fontsize = size_labels)
plt.yticks(Fontsize = size_labels)
for row in range(data['Model'].count()):
    if data['Species'][row] == 'Human':
        color = 'r'
    if data['Species'][row] == 'Dog':
        color = 'darkgreen'
    if data['Species'][row] == 'Rabbit':
        color = 'darkblue'
    if data['Species'][row] == 'GuineaPig':
        color = 'purple'
    if data['Species'][row] == 'Mammalian':
        color = 'orange'
    if data['Species'][row] == 'Rat':
        color = 'grey'
    if data['Species'][row] == 'Mouse':
        color = 'grey'


    if data['Tissue'][row] == 'Ventricle':
        shape = 'o'
    if data['Tissue'][row] == 'Atrium':
        shape = '+'
    if data['Tissue'][row] == 'Purkinje':
        shape = 'v'
    if data['Tissue'][row] == 'SAN':
        shape = 'x'  
     
    plt.scatter(data['Sodium_i_IC'][row], data['Potassium_i_IC'][row], color = color, marker = shape, s = 200)

"""
plt.subplot(1, 2, 2)
plt.title('At limit cycle', Fontsize = 30)
plt.xlabel('$[Na^+]_i$ (mM)', Fontsize = size_labels)
plt.ylabel('$[K^+]_i$ (mM)', Fontsize = size_labels)
plt.xlim([4, 19])
plt.ylim([50, 180])
plt.xticks([4, 6, 8, 10, 12, 14, 16, 18], [4, 6, 8, 10, 12, 14, 16, 18], Fontsize = size_labels)
plt.yticks(Fontsize = size_labels)
for row in range(data['Model'].count()-4):
    if data['Species'][row] == 'Human':
        color = 'r'
    if data['Species'][row] == 'Dog':
        color = 'darkgreen'
    if data['Species'][row] == 'Rabbit':
        color = 'darkblue'
    if data['Species'][row] == 'GuineaPig':
        color = 'purple'
    if data['Species'][row] == 'Mammalian':
        color = 'orange'
    if data['Species'][row] == 'Rat':
        color = 'grey'
    if data['Species'][row] == 'Mouse':
        color = 'grey'
        
    
    if data['Tissue'][row] == 'Ventricle':
        shape = 'o'
    if data['Tissue'][row] == 'Atrium':
        shape = '+'
    if data['Tissue'][row] == 'Purkinje':
        shape = 'v'
    if data['Tissue'][row] == 'SAN':
        shape = 'x'  
     
    plt.scatter(data['Sodium_i_steady'][row], data['Potassium_i_steady'][row], color = color, marker = shape, s = 200)

plt.tight_layout()

"""