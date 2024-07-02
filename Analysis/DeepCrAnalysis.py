#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Mar  5 13:10:59 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt
import glob

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

"""
AntennaPosition = 0 -1000 273500 ch0
AntennaPosition = 0 0 273500 ch1
AntennaPosition = 0 1000 273500 ch2
"""
z_ice = 2835 # m


Pos1 = np.array([0, -1000, 273500])
Pos2 = np.array([0, 0, 273500])
Pos3 = np.array([0, 1000, 273500])
AntPos = np.array([Pos1, Pos2, Pos3])/1e2
AntPos[:,2] -= z_ice

CoreasPath  = '/Users/chiche/Desktop/DeepCrAnalysis/Data/SIM000001/SIM000001_coreas'
GeantPath =  '/Users/chiche/Desktop/DeepCrAnalysis/Data/SIM000001/SIM000001_geant'


Ecgs_to_si = 3e10

CoreasFiles = glob.glob(CoreasPath + "/*")
GeantFiles = glob.glob(GeantPath + "/*")
Nant = len(CoreasFiles)

CoreasTraces = dict()
GeantTraces = dict()
EtotCoreas = np.zeros(Nant)
EtotGeant = np.zeros(Nant)

for i in range(Nant):
    
    CoreasTraces[i] = np.loadtxt(CoreasFiles[i], unpack = True)
    CoreasTraces[i][1:,:] *= Ecgs_to_si
    GeantTraces[i] = np.loadtxt(GeantFiles[i], unpack = True)
    GeantTraces[i][1:,:] *= Ecgs_to_si

# =============================================================================
#                       Antennas Positions
# =============================================================================

plt.scatter(AntPos[:, 1], AntPos[:, 2])
plt.xlabel("y [m]")
plt.ylabel("z [m]")
plt.ylim(-200, 0)
plt.title("Antennas positions (y,z) plane", fontsize =12)
plt.show()

plt.scatter(AntPos[:, 0], AntPos[:, 1])
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.title("Antennas positions (x,y) plane", fontsize =12)
plt.show()

# =============================================================================
#                           Efield Traces
# =============================================================================
# Channel k
k = 1

#Coreas
plt.plot(CoreasTraces[k][0,:]*1e9, CoreasTraces[k][1,:], label = "Ex")
plt.plot(CoreasTraces[k][0,:]*1e9, CoreasTraces[k][2,:], label = "Ey")
plt.plot(CoreasTraces[k][0,:]*1e9, CoreasTraces[k][3,:], label = "Ez")
plt.xlim(CoreasTraces[k][0,0]*1e9, 1200)
plt.legend()
plt.xlabel("Time [ns]")
plt.ylabel("E [$\mu V/m$]")
plt.title("CoREAS Traces Channel %.d" %k, fontsize =12)
plt.show()

for k in range(Nant):
    EtotCoreas[k]= np.sqrt(max(CoreasTraces[k][1,:]**2) + \
        max(CoreasTraces[k][2,:]**2) + max(CoreasTraces[k][3,:]**2))

#Geant
plt.plot(GeantTraces[k][0,:]*1e9, GeantTraces[k][1,:], label = "Ey")
plt.plot(GeantTraces[k][0,:]*1e9, GeantTraces[k][2,:], label = "Ey")
plt.plot(GeantTraces[k][0,:]*1e9, GeantTraces[k][3,:], label = "Ez")
plt.xlim(GeantTraces[k][0,0]*1e9, 560)
plt.legend()
plt.xlabel("Time [ns]")
plt.ylabel("E [$\mu V/m$]")
plt.title("Geant Traces Channel %.d" %k, fontsize =12)
plt.show()

for k in range(Nant):
    EtotGeant[k]= np.sqrt(max(GeantTraces[k][1,:]**2) + \
        max(GeantTraces[k][2,:]**2) + max(GeantTraces[k][3,:]**2))

plt.scatter(AntPos[:, 1], AntPos[:, 2], c= EtotCoreas, cmap ="jet")
plt.xlabel("y [m]")
plt.ylabel("z [m]")
plt.ylim(-200, 0)
cbar = plt.colorbar()
cbar.set_label("$E^{tot}_{coreas}$ [$\mu V/m$]")
plt.title("Antennas positions (y,z) plane", fontsize =12)
plt.show()

plt.scatter(AntPos[:, 1], AntPos[:, 2], c= EtotGeant, cmap ="jet")
plt.xlabel("y [m]")
plt.ylabel("z [m]")
cbar = plt.colorbar()
cbar.set_label("$E^{tot}_{geant}$ [$\mu V/m$]")
plt.title("Antennas positions (y,z) plane", fontsize =12)
plt.show()



plt.scatter(AntPos[:, 1], AntPos[:, 2], c= EtotCoreas, cmap ="jet")
plt.xlabel("y [m]")
plt.ylabel("z [m]")
plt.ylim(-200, 0)
cbar = plt.colorbar()
cbar.set_label("$E^{tot}_{coreas}$ [$\mu V/m$]")
plt.title("Antennas positions (y,z) plane", fontsize =12)
plt.show()

plt.scatter(AntPos[:, 1], AntPos[:, 2], c= EtotGeant/EtotCoreas, cmap ="jet")
plt.xlabel("y [m]")
plt.ylabel("z [m]")
cbar = plt.colorbar()
cbar.set_label("$E^{tot}_{geant}/E^{tot}_{coreas}$ [$\mu V/m$]")
plt.title("Antennas positions (y,z) plane", fontsize =12)
plt.show()


