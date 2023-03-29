# -*- coding: utf-8 -*-
"""
Created on Mon Mar 27 18:30:29 2023

@author: adria
"""

from matplotlib import pyplot as plt
import numpy as np

# Parámetros
# ========================================
datos=np.loadtxt('energy_values.txt') # Nombre del fichero de datos

t=datos[:,0]
p1=datos[:,1]
p2=datos[:,2]
p3=datos[:,3]
p4=datos[:,4]
p5=datos[:,5]
p6=datos[:,6]
p7=datos[:,7]
p8=datos[:,8]
p9=datos[:,9]

total_energy=[0.0]*len(t)

for i in range(len(t)):
    total_energy[i] = p1[i] + p2[i] + p3[i] + p4[i] + p5[i] + p6[i] + p7[i] + p8[i] + p9[i]

fig, ax = plt.subplots(1, 1, figsize=(40,10))

plt.rc('xtick', labelsize=16)
plt.rc('ytick', labelsize=16)

plt.subplot(1,1,1)
'''
plt.plot(t, p1, label='The Sun')
plt.plot(t, p2, label='Mercury')
plt.plot(t, p3, label='Venus')
plt.plot(t, p4, label='Earth')
plt.plot(t, p5, label='Mars')
plt.plot(t, p6, label='Jupiter')
plt.plot(t, p7, label='Saturn')
plt.plot(t, p8, label='Uranus')
plt.plot(t, p9, label='Neptune')
'''
plt.plot(t, total_energy, label='Total Energy ')
plt.xlabel('Time (days)', fontsize=20)
plt.ylabel('Energy', fontsize=20)
plt.title('Total energy of the system vs Time', fontsize=20)  
plt.show()
plt.legend(loc=1, prop={'size': 15})
plt.style.use('default')
