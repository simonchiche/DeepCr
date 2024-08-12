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
import pickle
from GetLayoutScaling import GetDepthScaling, GetEnergyScaling, Getradius
from scipy.interpolate import interp1d
import scipy
from FunctionsGetFluence import Norm, LoadTraces, GetPeakTraces, Traces_cgs_to_si, GetDepths, CorrectScaling, CombineTraces, CorrectLength, GetIntTraces, GetIntTracesSum, GetRadioExtent
from FunctionsPlotFluence import EfieldMap, PlotLDF, PlotTraces, plot_polarisation, PlotMaxTraces, PlotAllTraces, PlotLayer, PlotGivenTrace
from CleanCoreasTraces import CleanCoreasTraces
from scipy.interpolate import griddata
from datetime import datetime
from scipy.optimize import curve_fit


SimDir = "DeepCrLib"
SimName = "Rectangle_Proton_0.1_10_0_1"

date = datetime.today().strftime('%Y-%m-%d')
WorkPath = os.getcwd()
OutputPath = WorkPath + "/Plots/" + SimDir + "/" + date + "/" 
cmd = "mkdir -p " + OutputPath
p =subprocess.Popen(cmd, cwd=os.getcwd(), shell=True)
stdout, stderr = p.communicate()

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

Save = False

# =============================================================================
#                              Load Traces
# =============================================================================

Path =  "/Users/chiche/Desktop/DeepCrSearch"\
+ "/Simulations/" + SimDir + "/" + SimName + "/"
energy = float(Path.split("/")[-2].split("_")[2])
theta  = float(Path.split("/")[-2].split("_")[3])

SimDataPath = "/Users/chiche/Desktop/DeepCrSearch/Analysis/SimData/" + SimName 

if(not(os.path.exists(SimDataPath))):
    
    cmd = "mkdir -p " + SimDataPath
    p =subprocess.Popen(cmd, cwd=os.getcwd(), shell=True)
    stdout, stderr = p.communicate()
    Nant, Traces_C, Traces_G, Pos = LoadTraces(Path)
    
    np.save(SimDataPath + "/Nant", Nant)
    with open(SimDataPath + '/Traces_C.pkl', 'wb') as file:
        pickle.dump(Traces_C, file)
    with open(SimDataPath + '/Traces_G.pkl', 'wb') as file:
        pickle.dump(Traces_G, file)
    np.save(SimDataPath + "/Pos", Pos)
    
else:
    
    Nant  = np.load(SimDataPath + "/Nant.npy")
    with open(SimDataPath + '/Traces_C.pkl', 'rb') as file:
        Traces_C = pickle.load(file)
    with open(SimDataPath + '/Traces_G.pkl', 'rb') as file:
        Traces_G = pickle.load(file)
    Pos  = np.load(SimDataPath + "/Pos.npy", allow_pickle=True)
        
    
Nlay, Nplane, Depths = GetDepths(Pos)

# To resize the Coreas traces
CorrectLength(Traces_C, False)

#cgs to si
Traces_C = Traces_cgs_to_si(Traces_C)
Traces_G = Traces_cgs_to_si(Traces_G)


# =============================================================================
#                           Coherent sum
# =============================================================================

if(not(os.path.exists(SimDataPath + "/Traces_tot.pkl"))):
    Traces_tot = CombineTraces(Nant, Traces_C, Traces_G)
    with open(SimDataPath + '/Traces_tot.pkl', 'wb') as file:
        pickle.dump(Traces_C, file)   
else:
    with open(SimDataPath + '/Traces_tot.pkl', 'rb') as file:
        Traces_tot = pickle.load(file)    

# =============================================================================
#                           Get peak amplitude
# =============================================================================
# Peak value of the traces
ExC, EyC, EzC, EtotC = GetPeakTraces(Traces_C, Nant)
ExG, EyG, EzG, EtotG = GetPeakTraces(Traces_G, Nant)
Extot, Eytot, Eztot, Etot_peak = GetPeakTraces(Traces_tot, Nant)

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
PlotAllTraces(Nant, Traces_tot, 50, 1)

PlotGivenTrace(Traces_C, 310, "y")

# =============================================================================
#                         Compute Fluence
# =============================================================================

# correction if antenna with odd amplitude
if(max(Etot_int)>100*Etot_int[np.argsort(Etot_int)[1]]):
    kmax = np.argmax(Etot_int)
    Etot_int[kmax] = Etot_int[np.argsort(Etot_int)[1]]#np.mean(Etot_int)

EfieldMap(Pos, Nlay, Nplane, EtotC_int, "CoreasHilbert", False, energy, theta, OutputPath)

# Coreas Normalized
EfieldMap(Pos, Nlay, Nplane, EtotC_int/max(EtotC_int), "Coreas", False, energy, theta, OutputPath)

# Geant 
EfieldMap(Pos, Nlay, Nplane, EtotG_int, "GeantHilbert", True, energy, theta, OutputPath)
    
# Geant normalized
EfieldMap(Pos, Nlay, Nplane, EtotG_int/max(EtotG_int), "Geant", False, energy, theta, OutputPath)

# Total emission
EfieldMap(Pos, Nlay, Nplane, Etot_int, "Total", False, energy, theta, OutputPath)

#Total emission from peak
EfieldMap(Pos, Nlay, Nplane, np.maximum(EtotC, EtotG), "Total", False, energy, theta, OutputPath)

# Geant over CoREAS
EfieldMap(Pos, Nlay, Nplane, EtotG_int/EtotC_int, "GeantoverCoreas", False, energy, theta, OutputPath)
    

# =============================================================================
#                                 LDF
# =============================================================================

#Coreas
#PlotLDF(Pos, Nplane, EtotC, "Coreas", Nlay)
#Geant
#PlotLDF(Pos, Nplane, EtotG, "Geant", Nlay)

#PlotLayer(Pos, 2, Nplane, OutputPath)

# =============================================================================
#                            Polarisation
# =============================================================================
    
plot_polarisation(Pos[:Nplane,0], Pos[:Nplane,1], \
                  EtotC[:Nplane], -EyC[:Nplane], ExC[:Nplane], OutputPath)
    
# =============================================================================
#                    Parametrization of the spacing
# =============================================================================

radioextent, simextent, extent, maxpos, xminlay, xmaxlay = \
    GetRadioExtent(Nlay, Nplane, Pos, Etot_int)

print(radioextent/simextent)

### Plots

# radio extent and simulation extent vd depth
plt.plot(3216-np.array(Depths), radioextent, label = "radio")
plt.plot(3216-np.array(Depths), simextent, label = "sim")
plt.xlabel("Depth [m]")
plt.ylabel("Extent [m]")
plt.title("E =$%.2f\,$EeV, $\\theta=%.1f^{\circ}$" %(energy, theta), size =14)
plt.legend()
#plt.savefig(OutputPath + "LayoutExtent_E%.2f_th%.1f.pdf" \
 #            %(energy, theta), bbox_inches = "tight")
plt.show()

# filling factor
plt.scatter(3216-np.array(Depths), radioextent/simextent)
plt.xlabel("Depth [m]")
plt.ylabel("Filling factor [%]")
plt.title("E =$%.2f\,$EeV, $\\theta=%.1f^{\circ}$" %(energy, theta), size =14)
#plt.savefig(OutputPath + "FillingFactor_E%.2f_th%.1f.pdf" \
#             %(energy, theta), bbox_inches = "tight")
plt.show()



### Parametrization of the radio extent vs depth for different
# zenith angles and energies (to be saved in specific files) 
# see conclusions below
p0_0316 = np.array([31, 319, 336, 376, 407, 403, 461, 534, 538])
p1_0316 = np.array([2.813e-5, 4.1152e-4, 1.4380e-3, 1.9917e-4, 5.4562e-2,\
                    1.7811e-1, 2.5802e-1, 9.0471e-3, 6.1912e-2])
p2_0316 = np.array([3.07, 3.032, 2.77, 3.23, 1.932, 1.703, 1.664, 2.446, 2.0719])   

p0_01 = np.array([555, 553, 582, 645, 667, 789, 823, 922])
p1_01 = np.array([3.5956e-2, 8.0731e-2, 1.3556e-1, 9.4164e-2, 3.5645e-1, \
                  1.4066e-2, 3.7933e-1, 4.4415e-2])
p2_01 = np.array([2.04, 1.89, 1.842, 1.911, 1.691, 2.45, 1.7652, 2.27])   


r0 = Getradius(abs(np.array(Depths)-3216), p0_0316[0], p1_0316[0], p2_0316[0])
r1 = Getradius(abs(np.array(Depths)-3216), p0_0316[1], p1_0316[1], p2_0316[1])
r2 = Getradius(abs(np.array(Depths)-3216), p0_0316[2], p1_0316[2], p2_0316[2])
r3 = Getradius(abs(np.array(Depths)-3216), p0_0316[3], p1_0316[3], p2_0316[3])
r4 = Getradius(abs(np.array(Depths)-3216), p0_0316[4], p1_0316[4], p2_0316[4])
r5 = Getradius(abs(np.array(Depths)-3216), p0_0316[5], p1_0316[5], p2_0316[5])
r6 = Getradius(abs(np.array(Depths)-3216), p0_0316[6], p1_0316[6], p2_0316[6])
r7 = Getradius(abs(np.array(Depths)-3216), p0_0316[7], p1_0316[7], p2_0316[7])
r8 = Getradius(abs(np.array(Depths)-3216), p0_0316[8], p1_0316[8], p2_0316[8])

r0_17 = Getradius(abs(np.array(Depths)-3216), p0_01[0], p1_01[0], p2_01[0])
r1_17 = Getradius(abs(np.array(Depths)-3216), p0_01[1], p1_01[1], p2_01[1])
r2_17 = Getradius(abs(np.array(Depths)-3216), p0_01[2], p1_01[2], p2_01[2])
r3_17 = Getradius(abs(np.array(Depths)-3216), p0_01[3], p1_01[3], p2_01[3])
r4_17 = Getradius(abs(np.array(Depths)-3216), p0_01[4], p1_01[4], p2_01[4])
r5_17 = Getradius(abs(np.array(Depths)-3216), p0_01[5], p1_01[5], p2_01[5])
r6_17 = Getradius(abs(np.array(Depths)-3216), p0_01[6], p1_01[6], p2_01[6])
r7_17 = Getradius(abs(np.array(Depths)-3216), p0_01[7], p1_01[7], p2_01[7])


### Plots
# At 10^16.5 eV
plt.plot(abs(np.array(Depths)-3216), r1/p0_0316[1])
plt.plot(abs(np.array(Depths)-3216), r2/p0_0316[2])
plt.plot(abs(np.array(Depths)-3216), r3/p0_0316[3])
plt.plot(abs(np.array(Depths)-3216), r0/p0_0316[0])

plt.plot(abs(np.array(Depths)-3216), r4/p0_0316[4])
plt.plot(abs(np.array(Depths)-3216), r5/p0_0316[5])
plt.plot(abs(np.array(Depths)-3216), r6/p0_0316[6])
plt.plot(abs(np.array(Depths)-3216), r7/p0_0316[7])
plt.plot(abs(np.array(Depths)-3216), r8/p0_0316[8])
plt.xlabel("Depth [m]")
plt.ylabel("Scaling factor")
plt.title("E = $10^{16.5}$ eV", fontsize =14)
#plt.savefig(OutputPath + "Scaling_vs_Depth_E%.2f.pdf" \
#             %(energy), bbox_inches = "tight")
plt.show()

# At 10^17 eV
plt.plot(abs(np.array(Depths)-3216), r0_17/p0_01[0])
plt.plot(abs(np.array(Depths)-3216), r1_17/p0_01[1])
plt.plot(abs(np.array(Depths)-3216), r2_17/p0_01[2])
plt.plot(abs(np.array(Depths)-3216), r3_17/p0_01[3])
plt.plot(abs(np.array(Depths)-3216), r4_17/p0_01[4])
plt.plot(abs(np.array(Depths)-3216), r5_17/p0_01[5])
plt.plot(abs(np.array(Depths)-3216), r6_17/p0_01[6])
plt.plot(abs(np.array(Depths)-3216), r7_17/p0_01[7])
plt.xlabel("Depth [m]")
plt.ylabel("Scaling factor")
plt.title("E = $10^{17}$ eV", fontsize =14)
plt.savefig(OutputPath + "Scaling_vs_Depth_E%.2f.pdf" \
             %(energy), bbox_inches = "tight")
plt.show()

# 10^16.5 eV and 10^17 eV mixed

plt.plot(abs(np.array(Depths)-3216), r1/p0_0316[1], linestyle =":")
plt.plot(abs(np.array(Depths)-3216), r2/p0_0316[2], linestyle =":")
plt.plot(abs(np.array(Depths)-3216), r3/p0_0316[3], linestyle =":")
plt.plot(abs(np.array(Depths)-3216), r0/p0_0316[0], linestyle =":")

plt.plot(abs(np.array(Depths)-3216), r4/p0_0316[4], linestyle =":")
plt.plot(abs(np.array(Depths)-3216), r5/p0_0316[5], linestyle =":")
plt.plot(abs(np.array(Depths)-3216), r6/p0_0316[6], linestyle =":")
plt.plot(abs(np.array(Depths)-3216), r7/p0_0316[7], linestyle =":")
plt.plot(abs(np.array(Depths)-3216), r8/p0_0316[8], linestyle =":")

plt.plot(abs(np.array(Depths)-3216), r0_17/p0_01[0])

plt.plot(abs(np.array(Depths)-3216), r1_17/p0_01[1])
plt.plot(abs(np.array(Depths)-3216), r2_17/p0_01[2])
plt.plot(abs(np.array(Depths)-3216), r3_17/p0_01[3])
plt.plot(abs(np.array(Depths)-3216), r4_17/p0_01[4])
plt.plot(abs(np.array(Depths)-3216), r5_17/p0_01[5])
plt.plot(abs(np.array(Depths)-3216), r6_17/p0_01[6])
plt.plot(abs(np.array(Depths)-3216), r7_17/p0_01[7])
plt.show()


# CCL: the seventh parametrization (in index) at 10^16.5 eV seems to reproduce well the dependency 
# of the radio extent with the depth independently of the simulation energy and depth

"""
# preliminary test for analytic modeling

theta16 = np.array([0, 10, 20, 28, 34, 39, 43, 47, 50])
theta17 = np.array([10, 20, 28, 34, 39, 43, 47, 50])

cos_theta16 = np.cos(theta16*np.pi/180.0)
cos_theta17 = np.cos(theta17*np.pi/180.0)
scalingAll = np.concatenate((cos_theta16, cos_theta17))
extentmax =np.array([r0[0]/p0_0316[0], r1[0]/p0_0316[1], r2[0]/p0_0316[2], \
   r3[0]/p0_0316[3], r4[0]/p0_0316[4], r5[0]/p0_0316[5], r6[0]/p0_0316[6], \
   r7[0]/p0_0316[7], r8[0]/p0_0316[8], r0_17[0]/p0_01[0], r1_17[0]/p0_01[1], r2_17[0]/p0_01[2], r3_17[0]/p0_01[3], r4_17[0]/p0_01[4], r5_17[0]/p0_01[5], r6_17[0]/p0_01[6], r7_17[0]/p0_01[7]])
plt.scatter(scalingAll, extentmax)

popt, pcov = curve_fit(Getradius, abs(np.array(Depths)-3216), radioextent)

r = Getradius(abs(np.array(Depths)-3216), popt[0], popt[1], popt[2])
rtest = Getradius(abs(np.array(Depths)-3216), p0_0316[2], p1_0316[2], p2_0316[2])

plt.scatter(abs(np.array(Depths)-3216), radioextent/popt[0])
plt.plot(abs(np.array(Depths)-3216), rtest/p0_0316[2])
plt.plot(abs(np.array(Depths)-3216), r/popt[0])
plt.show()
"""

