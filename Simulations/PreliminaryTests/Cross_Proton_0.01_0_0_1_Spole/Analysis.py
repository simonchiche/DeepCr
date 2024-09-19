#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 17 18:37:55 2024

@author: chiche
"""
#3216

import numpy as np
import glob
import matplotlib.pyplot as plt


CoreasFile = glob.glob("/Users/chiche/Desktop/Cross_Proton_0.01_0_0_1/Coreas/raw*")
GeantFile = glob.glob("/Users/chiche/Desktop/Cross_Proton_0.01_0_0_1/Geant/antenna*")
PosFile = glob.glob("/Users/chiche/Desktop/Cross_Proton_0.01_0_0_1/Parameters/*.list")


# Initialize empty lists for x, y, and z coordinates
x = []
y = []
z = []

# Open the file
with open(PosFile[0], "r") as file:
    # Iterate through each line in the file
    for line in file:
        # Split the line by whitespace
        parts = line.split()
        # Extract the x, y, and z coordinates from the line and convert them to floats
        x_coord = float(parts[2])
        y_coord = float(parts[3])
        z_coord = float(parts[4])
        # Append the coordinates to their respective lists
        x.append(x_coord)
        y.append(y_coord)
        z.append(z_coord)

x = np.array(x)
y = np.array(y)
z = np.array(z)


plt.scatter(x,y)
plt.show()

CoreasE = []
GeantE = []


for i in range(len(CoreasFile)):
    CoreasE.append(np.loadtxt(CoreasFile[i], unpack = True))
    GeantE.append(np.loadtxt(GeantFile[i], unpack = True))

CoreasE = np.array(CoreasE)
GeantE = np.array(GeantE)


EtotC = np.zeros(len(CoreasE))
EtotG = np.zeros(len(GeantE))

for i in range(len(EtotC)):
    
    EtotC[i]= max(np.sqrt(CoreasE[i][1,:]**2 + CoreasE[i][2,:]**2 + CoreasE[i][3,:]**2))
    
    EtotG[i]= max(np.sqrt(GeantE[i][1,:]**2 + GeantE[i][2,:]**2 + GeantE[i][3,:]**2))
    
    


plt.scatter(x,y, c = EtotC)

plt.plot(CoreasE[10][0,:], CoreasE[10][2,:])
plt.xlim(0.4e-6,0.8e-6)