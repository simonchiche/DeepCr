#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 02:21:35 2024

@author: chiche
"""


import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt
import glob
import sys
from GetLayoutScaling import GetDepthScaling, GetEnergyScaling
from scipy.interpolate import interp1d
import scipy
from FunctionsGetFluence import Norm, LoadTraces, GetPeakTraces, Traces_cgs_to_si, GetDepths, CorrectScaling, CombineTraces, CorrectLength, GetIntTraces, GetIntTracesSum
from FunctionsPlotFluence import EfieldMap, PlotLDF, PlotTraces, plot_polarisation, PlotMaxTraces, PlotAllTraces, PlotLayer
from CleanCoreasTraces import CleanCoreasTraces
from scipy.interpolate import griddata
from datetime import datetime

date = datetime.today().strftime('%Y-%m-%d')
WorkPath = os.getcwd()
OutputPath = WorkPath + "/Plots/" + date + "/"
cmd = "mkdir -p " + OutputPath
p =subprocess.Popen(cmd, cwd=os.getcwd(), shell=True)
stdout, stderr = p.communicate()

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

# =============================================================================
#                              Load Traces
# =============================================================================

Path =  "/Users/chiche/Desktop/DeepCrSearch"\
+ "/Simulations/SpoleDoubleRing/Rectangle_Proton_0.0316_0_0_1/"
energy = float(Path.split("/")[-2].split("_")[2])
theta  = float(Path.split("/")[-2].split("_")[3])

Nant, Traces_C, Traces_G, Pos = LoadTraces(Path)
Nlay, Nplane, Depths = GetDepths(Pos)

# To resize the Coreas traces
CorrectLength(Traces_C, False)

#cgs to si
Traces_C = Traces_cgs_to_si(Traces_C)
Traces_G = Traces_cgs_to_si(Traces_G)


# =============================================================================
#                           Coherent sum
# =============================================================================

Traces_tot = CombineTraces(Nant, Traces_C, Traces_G)

# =============================================================================
#                           Get peak amplitude
# =============================================================================
# Peak value of the traces
ExC, EyC, EzC, EtotC = GetPeakTraces(Traces_C, Nant)
ExG, EyG, EzG, EtotG = GetPeakTraces(Traces_G, Nant)
Extot, Eytot, Eztot, Etot_peak = GetPeakTraces(Traces_tot, Nant)

#np.savetxt("EfieldCoreas_E%.1f_th%.1f.txt" %(energy, theta), np.array([ExC, EyC, EzC]))
#np.savetxt("EfieldGeant_E%.1f_th%.1f.txt" %(energy, theta), np.array([ExG, EyG, EzG]))


# =============================================================================
#                                 Get integral
# =============================================================================

ExC_int, EyC_int, EzC_int, EtotC_int = GetIntTraces(Traces_C, Nant)

ExG_int, EyG_int, EzG_int, EtotG_int = GetIntTraces(Traces_G, Nant)

Ex_tot_int, Ey_tot_int, Ez_tot_int, Etot_int = GetIntTraces(Traces_tot, Nant)

# =============================================================================
#                             Plot Traces
# =============================================================================

### Plot max traces
# Geant
PlotMaxTraces(Traces_G, EtotG, 1)
#Coreas
PlotMaxTraces(Traces_C, EtotC, 1)

# Plot all traces above a given threshold
PlotAllTraces(Nant, Traces_tot, 50, 10)

i = 310
plt.plot(Traces_tot[i][:, 0]*1e9,Traces_tot[i][:, 2])
plt.xlabel("Time [ns]")
plt.ylabel("E [$\mu V/m$]")
plt.title(r"Proton  E = 0.01 EeV - $\theta = 0^{\circ}$ $\phi = 0^{\circ}$")
#plt.savefig(OutputPath + "Double_Pulse1.pdf", bbox_inches = "tight")
plt.show()


plt.plot(Traces_C[4*Nplane +337][:, 0]*1e9, -Traces_C[4*Nplane + 337][:, 2], label = "Depth = 0 m")
plt.plot(Traces_tot[337][:, 0]*1e9,Traces_tot[337][:, 2], label = "Depth = 100 m")
plt.xlabel("Time [ns]")
plt.ylabel("E [$\mu V/m$]")
plt.title(r"Proton  E = 0.032 EeV - $\theta = 0^{\circ}$ $\phi = 0^{\circ}$")
plt.legend()
#plt.savefig("/Users/chiche/Desktop/Tpulse_VS_depth.pdf", bbox_inches = "tight")
plt.show()

for j in range(337, 338, 1):
    print(j)
    plt.plot(Traces_C[4*Nplane + j][:, 0]*1e9,Traces_C[4*Nplane + j][:, 2])
    plt.xlabel("Time [ns]")
    plt.ylabel("E [$\mu V/m$]")
    plt.title(r"Proton  E = 0.01 EeV - $\theta = 0^{\circ}$ $\phi = 0^{\circ}$")
    #plt.savefig(OutputPath + "Double_Pulse1.pdf", bbox_inches = "tight")
    plt.show()

# =============================================================================
#                         Compute Fluence
# =============================================================================

### if antenna in (0,0) we remove it
#EtotG[3280]= 100
#r = np.sqrt(Pos[:,0]**2 + Pos[:,1]**2)
#EtotG[r<100]=0


# Coreas
EfieldMap(Pos, Nlay, Nplane, EtotC, "CoreasPeakSpole", True, energy, theta, OutputPath)

# Coreas Normalized
EfieldMap(Pos, Nlay, Nplane, EtotC_int/max(EtotC_int), "Coreas", False, energy, theta, OutputPath)

# Geant 
EfieldMap(Pos, Nlay, Nplane, EtotG_int, "GeantHilbertSpole", True, energy, theta, OutputPath)
    
# Geant normalized
EfieldMap(Pos, Nlay, Nplane, EtotG_int/max(EtotG_int), "Geant", False, energy, theta, OutputPath)

# Total emission
EfieldMap(Pos, Nlay, Nplane, Etot_int, "Total", True, energy, theta, OutputPath)

#Total emission from peak
EfieldMap(Pos, Nlay, Nplane, np.maximum(EtotC, EtotG), "Total", False, energy, theta, OutputPath)



# =============================================================================
#                                 LDF
# =============================================================================

#Coreas
PlotLDF(Pos, Nplane, EtotC, "Coreas", Nlay)
#Geant
PlotLDF(Pos, Nplane, EtotG, "Geant", Nlay)

k = 2
plt.plot(Pos[:int(Nplane/2),0], EtotG[k*Nplane:(k*Nplane+int(Nplane/2))])
plt.show()

PlotLayer(Pos, k, Nplane, OutputPath)

# =============================================================================
#                            Polarisation
# =============================================================================
    
plot_polarisation(Pos[:Nplane,0], Pos[:Nplane,1], \
                  EtotC[:Nplane], -EyC[:Nplane], ExC[:Nplane], OutputPath)
    
# =============================================================================
#                          interpolation
# =============================================================================
"""

grid_x, grid_y = \
np.mgrid[min(Pos[:Nplane,0]):max(Pos[:Nplane,1]):100j, min(Pos[:Nplane,0]):max(Pos[:Nplane,1]):100j]
#grid_x, grid_y = np.mgrid[-1:1:300j, -1:1:300j]
points = np.array([Pos[:Nplane,0], Pos[:Nplane,1]]).T

# On associe les valeurs du champ interpolé à chaque point de la grille
grid_z1 = griddata(points, EtotG_int[:Nplane], (grid_x, grid_y), method='cubic')


plt.pcolormesh(grid_x, grid_y, grid_z1, cmap ="jet")
cbar = plt.colorbar()
plt.xlim(-250,450)
plt.ylim(-220,350)
plt.show()

"""
