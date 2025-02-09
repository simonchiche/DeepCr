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
sys.path.append("/Users/chiche/Desktop/DeepCrSearch/Analysis/")
from Modules.SimParam.GetLayoutScaling import GetDepthScaling, GetEnergyScaling
from scipy.interpolate import interp1d
import scipy
from Modules.Fluence.FunctionsGetFluence import  Norm, LoadTraces, GetPeakTraces, Traces_cgs_to_si, GetDepths, CorrectScaling, CombineTraces, CorrectLength, GetIntTraces, GetIntTracesSum, GetRadioExtent
from Modules.Fluence.FunctionsPlotFluence import EfieldMap, PlotLDF, PlotTraces, plot_polarisation, PlotMaxTraces, PlotAllTraces, PlotLayer, PlotGivenTrace, PlotAllChannels
from CleanCoreasTraces import CleanCoreasTraces
from Modules.SimParam.PlotRadioSimExtent import PlotFillingFactor, PlotRadioSimExtent
from scipy.interpolate import griddata
from datetime import datetime
from scipy.optimize import curve_fit
from scipy.signal import butter, filtfilt
from Modules.Fluence.FunctionsRadiationEnergy import GetFluence, GetRadiationEnergy

PowerDataPath = "/Users/chiche/Desktop/DeepCrSearch/Analysis/GetFluence/Data/Power/"

SimDir =     "DeepCrLib" #"InterpSim"
SimName = "Rectangle_Proton_0.316_50_0_1"

# We create a directory where the outputs are stored
date = datetime.today().strftime('%Y-%m-%d')
WorkPath = os.getcwd()
OutputPath = WorkPath + "/Plots/" + SimDir + "/" + date + "/" 
#print(OutputPath)
#sys.exit()
cmd = "mkdir -p " + OutputPath
p =subprocess.Popen(cmd, cwd=os.getcwd(), shell=True)
stdout, stderr = p.communicate()

font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)
zenith = float(SimName.split("_")[3])
Save = False

# =============================================================================
#                              Load Traces
# =============================================================================

Path =  "/Users/chiche/Desktop/DeepCrSearch"\
+ "/Simulations/" + SimDir + "/" + SimName + "/"
energy = float(Path.split("/")[-2].split("_")[2])
theta  = float(Path.split("/")[-2].split("_")[3])

SimDataPath = "/Users/chiche/Desktop/DeepCrSearch/Analysis/SimData/" + SimDir + "/" + SimName 

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

# To resize the length of the Coreas traces
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
#                                 Write Power
# =============================================================================

spacing = abs(Pos[1,0] - Pos[0,0])
Eair = np.sum(EtotC_int[Pos[:,2] == 3116])
Eice = np.sum(EtotG_int[Pos[:,2] == 3116])
'''
fair = np.sum(EtotC_int[Pos[:,2] == 3116]**2)#*spacing**2
fice = np.sum(EtotG_int[Pos[:,2] == 3116]**2)#*spacing**2

with open("fluencevs_theta.txt", "a") as f:
    f.write(f"{zenith}\t{fair}\t{fice}\n")
'''
with open(PowerDataPath + "Efield_vs_theta.txt", "a") as f:
    f.write(f"{zenith}\t{Eair}\t{Eice}\n")

sys.exit()

# =============================================================================
#                             Plot Traces
# =============================================================================

### Plot max traces
# Geant
##PlotMaxTraces(Traces_G, EtotG, 1)
#Coreas
##PlotMaxTraces(Traces_C, EtotC, 5)

# Plot all traces above a given threshold
##PlotAllTraces(Nant, Traces_tot, 100, 5)

##PlotGivenTrace(Traces_C, 310, "y")

##PlotAllChannels(Traces_C, 1068)
##PlotAllChannels(Traces_G, 1068)

# =============================================================================
#                         Compute Fluence
# =============================================================================

# correction if antenna with odd amplitude
if(max(Etot_int)>100*Etot_int[np.argsort(Etot_int)[1]]):
    kmax = np.argmax(Etot_int)
    Etot_int[kmax] = Etot_int[np.argsort(Etot_int)[1]]#np.mean(Etot_int)


EfieldMap(Pos, Nlay, Nplane, EtotC_int, "CoreasHilbert", \
          True, energy, theta, OutputPath)

# Coreas Normalized
EfieldMap(Pos, Nlay, Nplane, EtotC_int/max(EtotC_int), "Coreas",\
          False, energy, theta, OutputPath)

# Geant 
EfieldMap(Pos, Nlay, Nplane, EtotG_int, "GeantHilbert",\
          True, energy, theta, OutputPath)
    
# Geant normalized
EfieldMap(Pos, Nlay, Nplane, EtotG_int/max(EtotG_int), "Geant", \
          False, energy, theta, OutputPath)

# Total emission
EfieldMap(Pos, Nlay, Nplane, Etot_int, "Total", \
          False, energy, theta, OutputPath)

#Total emission from peak
#EfieldMap(Pos, Nlay, Nplane, np.maximum(EtotC, EtotG), "Total",\
#          False, energy, theta, OutputPath)

# Geant over CoREAS
EfieldMap(Pos, Nlay, Nplane, EtotG_int/EtotC_int, "GeantoverCoreas",\
          False, energy, theta, OutputPath)
    

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
                  EtotC_int[:Nplane], -EyC_int[:Nplane], ExC_int[:Nplane], OutputPath)
    
# =============================================================================
#                    Parametrization of the spacing
# =============================================================================

radioextent, simextent, extent, maxpos, xminlay, xmaxlay = \
    GetRadioExtent(Nlay, Nplane, Pos, Etot_int)

print(radioextent/simextent)

## Plots

PlotRadioSimExtent(Depths, radioextent, simextent)
PlotFillingFactor(Depths, radioextent, simextent)



def bandpass_filter(signal, fs, lowcut, highcut, order=4):
    """
    Apply a bandpass filter to a signal.
    
    Parameters:
    - signal: array-like, the input signal (E(t)).
    - fs: float, the sampling frequency of the signal in Hz.
    - lowcut: float, the lower bound of the frequency band in Hz.
    - highcut: float, the upper bound of the frequency band in Hz.
    - order: int, the order of the Butterworth filter (default is 4).
    
    Returns:
    - filtered_signal: array-like, the filtered signal.
    """
    nyquist = 0.5 * fs  # Nyquist frequency
    low = lowcut / nyquist
    high = highcut / nyquist

    # Design the Butterworth bandpass filter
    b, a = butter(order, [low, high], btype='band')

    # Apply the filter using filtfilt for zero phase shift
    filtered_signal = filtfilt(b, a, signal)
    return filtered_signal


fs = 5e9
lowcut = 50e6  # Lower bound of the frequency band in Hz (50 MHz)
highcut = 2e9  # Upper bound of the frequency band in Hz (2000 MHz)

Eg_f = np.zeros(len(Traces_G))
filtered_signal = []
for i in range(len(Traces_G)):
    Exg_f = bandpass_filter(Traces_G[i][:,1], fs, lowcut, highcut)
    Eyg_f = bandpass_filter(Traces_G[i][:,2], fs, lowcut, highcut)
    Ezg_f = bandpass_filter(Traces_G[i][:,3], fs, lowcut, highcut)
    
    Etotg_f = np.sqrt(Exg_f**2 + Eyg_f**2 + Ezg_f**2)
    #filtered_signal.append(bandpass_filter(Traces_G[i][:,2], fs, lowcut, highcut))
    Eg_f[i] = max(abs(Etotg_f))
    
#filtered_signal =np.array(filtered_signal)



# Geant 
EfieldMap(Pos, Nlay, Nplane, Eg_f, "GeantHilbert",\
          True, energy, theta, OutputPath)





