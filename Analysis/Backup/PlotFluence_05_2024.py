#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 02:21:35 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt
import glob
import sys
from GetLayoutScaling import GetDepthScaling, GetEnergyScaling
from scipy.interpolate import interp1d

from FunctionsGetFluence import Norm, LoadTraces, GetPeakTraces, Traces_cgs_to_si, GetDepths, CorrectScaling
from FunctionsPlotFluence import EfieldMap, PlotLDF, PlotTraces
from CleanCoreasTraces import CleanCoreasTraces

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

# =============================================================================
#                              Load Traces
# =============================================================================

Path =  "/Users/chiche/Desktop/DeepCrSearch"\
+ "/Simulations/DeepCrLib/Rectangle_Proton_0.0316_10_0_1/"
energy = float(Path.split("/")[-2].split("_")[2])
theta  = float(Path.split("/")[-2].split("_")[3])

Nant, Traces_C, Traces_G, Pos = LoadTraces(Path)
Nlay, Nplane, Depths = GetDepths(Pos)

CorrectLength = False
if(CorrectLength):
    for i in range(len(Traces_C)):
        if(i%10==0): print(i)
        np.savetxt("/Users/chiche/Desktop/DeepCrSearch"\
         + "/Simulations/DeepCrLib/Rectangle_Proton_0.0316_10_0_1/NewCoreas/raw_ch%.d.dat" %i, Traces_C[i][:5000,:])



#cgs to si
#Traces_C = Traces_cgs_to_si(Traces_C)
Traces_G = Traces_cgs_to_si(Traces_G)

"""
def align_and_sum_signals(t1, E1, t2, E2, method='linear'):
    # Determine the overall start and end time
    start_time = min(t1[0], t2[0])
    end_time = max(t1[-1], t2[-1])
    
    # Determine the time bin size (assuming it is the same for both)
    time_bin = np.diff(t1).mean()  # assuming t1 and t2 have the same time bin
    aligned_time = np.arange(start_time, end_time + time_bin, time_bin)
    
    # Interpolate E1 and E2 to the aligned_time grid
    interp_E1 = interp1d(t1, E1, kind=method, fill_value=0, bounds_error=False)
    interp_E2 = interp1d(t2, E2, kind=method, fill_value=0, bounds_error=False)
    
    # Evaluate the interpolated signals on the common time grid
    aligned_E1 = interp_E1(aligned_time)
    aligned_E2 = interp_E2(aligned_time)
    
    # Sum the interpolated signals
    ##summed_signal = aligned_E1 + aligned_E2
    
    return aligned_time, aligned_E1, aligned_E2  #summed_signal

Traces_C_full = dict()
Traces_G_full = dict()
TracesAll = dict()

for i in range(len(Traces_C)):
    if(i%10==0): print(i)
    
    t, Ex_c, Ex_g = align_and_sum_signals(Traces_C[i][:,0], Traces_C[i][:,1],\
                                          Traces_G[i][:,0], Traces_G[i][:,1], method='linear')
    t, Ey_c, Ey_g = align_and_sum_signals(Traces_C[i][:,0], Traces_C[i][:,2],\
                                          Traces_G[i][:,0], Traces_G[i][:,2], method='linear')
    t, Ez_c, Ez_g = align_and_sum_signals(Traces_C[i][:,0], Traces_C[i][:,3],\
                                          Traces_G[i][:,0], Traces_G[i][:,3], method='linear')
    
    ExAll = Ex_c + Ex_g
    EyAll = Ey_c + Ey_g
    EzAll = Ez_c + Ez_g
    
    Traces_C_full[i] = np.array([t, Ex_c, Ey_c, Ez_c])
    Traces_G_full[i] = np.array([t, Ex_g, Ey_g, Ez_g])
    TracesAll[i] = np.array([t, ExAll, EyAll, EzAll])

"""
# =============================================================================
#                                 Plot Traces
# =============================================================================

# We reduce the size of the Coreas traces (long track with no signal)
#Traces_C = CleanCoreasTraces(Traces_C)

## We plot the traces
## PlotTraces(Traces_G, 0, 1)


### We plot the highest energy traces in priority

# Peak value of the traces
ExC, EyC, EzC, EtotC = GetPeakTraces(Traces_C, Nant)
ExG, EyG, EzG, EtotG = GetPeakTraces(Traces_G, Nant)
#Exall, Eyall, Ezall, Etotall = GetPeakTraces(TracesAll, Nant)


### Geant
# Plot max traces
MaxId = np.argsort(abs(EtotG))
NantPlot = 10
for i in range(NantPlot):
    arg = MaxId[-(i+1)]
    PlotTraces(Traces_G, arg, arg+1)

# Coreas
MaxId = np.argsort(abs(EtotC))#np.arange(0,100,1)#
NantPlot = 10
for i in range(NantPlot):
    arg = MaxId[-(i+1)]
    PlotTraces(Traces_C, arg, arg+1)
    

PlotAll = False  
NantPlot = Nant  
if(PlotAll):
    for i in range(NantPlot):
        #arg = MaxId[-(i+1)]
        if(max(abs(Traces_C[i][:, 1]))>50):
            PlotTraces(Traces_C, i, i+1)

"""
EtotAll = np.zeros(len(EtotC))
for i in range(len(EtotC)):
    EtotAll[i] = max(EtotC[i], EtotG[i])
"""  
# =============================================================================
#                         Compute Fluence
# =============================================================================

### if antenna in (0,0) we remove it
EtotG[3280]= 100
# Coreas
EfieldMap(Pos, Nlay, Nplane, EtotC, "Coreas", False, energy, theta)

r = np.sqrt(Pos[:,0]**2 + Pos[:,1]**2)
EtotG[r<100]=0
EfieldMap(Pos, Nlay, Nplane, np.maximum(EtotC, EtotG), "Total", True, energy, theta)


# Coreas Normalized
EfieldMap(Pos, Nlay, Nplane, EtotC/max(EtotC), "Coreas", False, energy, theta)

# Geant 
EfieldMap(Pos, Nlay, Nplane, EtotG, "Geant", False, energy, theta)
    
# Geant normalized
EfieldMap(Pos, Nlay, Nplane, EtotG/max(EtotG), "Geant", False, energy, theta)

# =============================================================================
#                                 LDF
# =============================================================================

#Coreas
PlotLDF(Pos, Nplane, EtotC, "Coreas", Nlay)
#Geant
PlotLDF(Pos, Nplane, EtotG, "Geant", Nlay)

#180 antennes par dag file


k = 2
plt.scatter(Pos[:int(Nplane/2),0], EtotG[k*Nplane:(k*Nplane+int(Nplane/2))])
plt.show()

plt.scatter(Pos[4*Nplane:5*Nplane,0], Pos[4*Nplane:5*Nplane,1])
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.tight_layout()
#plt.savefig("/Users/chiche/Desktop/RectGrid.pdf")
plt.show()



def plot_polarisation(vxb, vxvxb, Etot_sp, Evxb, Evxvxb):
    
    # function that plots the normalised polarisation in a given plane
    
    #r = np.sqrt(Evxb**2 + Evxvxb**2)
    plt.scatter(vxb, vxvxb, color = "white")
    #cbar = plt.colorbar()
    plt.xlabel('k x b [m]')
    plt.ylabel('k x k x b [m]')
    #plt.xlim(-200,200)
    #plt.ylim(-200,200)
    #plt.title(r'$\phi$ $= %.2f \degree$, $\theta$ $= %.2f \degree$, E = %.3f Eev' %(azimuth, zenith, energy), fontsize = 12)
    #cbar.set_label(r"$ E\ [\mu V/m]$")
    plt.quiver(vxb, vxvxb, Evxb/2000, Evxvxb/2000, scale =1, width=0.005)
    plt.xlim(-120,120)
    plt.ylim(-120,120)
    plt.tight_layout()
    plt.savefig('/Users/chiche/Desktop/polarisation_sp.png', dpi = 500)
    plt.show()
    

    
plot_polarisation(Pos[:Nplane,0], Pos[:Nplane,1], EtotC[:Nplane], -EyC[:Nplane], ExC[:Nplane])
    
