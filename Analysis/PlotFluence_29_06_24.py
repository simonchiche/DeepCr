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
from GetLayoutScaling import GetDepthScaling, GetEnergyScaling
from scipy.interpolate import interp1d
import scipy
from FunctionsGetFluence import Norm, LoadTraces, GetPeakTraces, Traces_cgs_to_si, GetDepths, CorrectScaling, CombineTraces, CorrectLength, GetIntTraces, GetIntTracesSum
from FunctionsPlotFluence import EfieldMap, PlotLDF, PlotTraces, plot_polarisation, PlotMaxTraces, PlotAllTraces, PlotLayer, PlotGivenTrace
from CleanCoreasTraces import CleanCoreasTraces
from scipy.interpolate import griddata
from datetime import datetime

SimDir = "DeepCrLib"
SimName = "Rectangle_Proton_0.1_39_0_1"

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
    with open(SimDataPath + '/Traces_G.pkl', 'rb') as file:
        Traces_G = pickle.load(file)    

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

# Coreas
##EtotC_int[364] = 0
##EtotC_int[1093] = 0

EfieldMap(Pos, Nlay, Nplane, EtotC_int, "CoreasHilbert", True, energy, theta, OutputPath)

# Coreas Normalized
EfieldMap(Pos, Nlay, Nplane, EtotC_int/max(EtotC_int), "Coreas", True, energy, theta, OutputPath)

# Geant 
EfieldMap(Pos, Nlay, Nplane, EtotG_int, "GeantHilbert", True, energy, theta, OutputPath)
    
# Geant normalized
EfieldMap(Pos, Nlay, Nplane, EtotG_int/max(EtotG_int), "Geant", True, energy, theta, OutputPath)

# Total emission
EfieldMap(Pos, Nlay, Nplane, Etot_int, "Total", Save, energy, True, OutputPath)

#Total emission from peak
EfieldMap(Pos, Nlay, Nplane, np.maximum(EtotC, EtotG), "Total", Save, energy, theta, OutputPath)

# Geant 
EfieldMap(Pos, Nlay, Nplane, EtotG_int/EtotC_int, "GeantoverCoreas", Save, energy, theta, OutputPath)
    

# =============================================================================
#                                 LDF
# =============================================================================

#Coreas
PlotLDF(Pos, Nplane, EtotC, "Coreas", Nlay)
#Geant
PlotLDF(Pos, Nplane, EtotG, "Geant", Nlay)

PlotLayer(Pos, 2, Nplane, OutputPath)

# =============================================================================
#                            Polarisation
# =============================================================================
    
plot_polarisation(Pos[:Nplane,0], Pos[:Nplane,1], \
                  EtotC[:Nplane], -EyC[:Nplane], ExC[:Nplane], OutputPath)
    

