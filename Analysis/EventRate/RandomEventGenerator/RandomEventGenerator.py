#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 14:38:32 2024

@author: chiche
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from GetRandomDirection import GetRandomTheta, GetRandomPhi
from GenRandomE import GetRandomE
from CreateRandomEvent import GetDeepAntennaLayerEvents, SaveEvent
from RunHDF5converter import RunHDF5converter

#SimulationPath
SimPath = "/Users/chiche/Desktop/DeepCrSearch/Simulations/DeepCrLib"
DataPath = "/Users/chiche/Desktop/DeepCrSearch/Analysis/SimData"
DataFilesPath = "/Users/chiche/Desktop/REGevents/SimFiles/"
LibPath = "/Users/chiche/Desktop/REGevents/Events"
HDF5path = "/Users/chiche/Desktop/REGevents/HDF5files"

# ground level in meters
glevel =3216

# Number of events 
Nev = 1000

# zenith limits
min_th = 0
max_th = 50
theta_rand = GetRandomTheta(min_th, max_th, Nev)

# azimuth limits 
min_phi = 0
max_phi = 360  
phi_rand = GetRandomPhi(min_phi, max_phi, Nev)

# Energy limits (log)
pmin = 16
pmax = 17.5
Erand = GetRandomE(pmin, pmax, Nev)

for i in range(len(Erand)):
    
    #AntPow, Traces_C_pow, Traces_G_pow, xrand, yrand = \
    #SelectPowerLineEvent(Erand[i], theta_rand[i], phi_rand[i], glevel, SimPath, DataPath)
    
    AntPow, Traces_C_pow, Traces_G_pow, xrand, yrand, Selfile = \
    GetDeepAntennaLayerEvents(0.1, 10, 24, glevel, SimPath, DataPath)
        
    savedir = \
    SaveEvent(AntPow, Traces_C_pow, Traces_G_pow, Selfile, LibPath, DataFilesPath, i)
    
    print("saved")
    
    RunHDF5converter(HDF5path, savedir)

    sys.exit()
   
        