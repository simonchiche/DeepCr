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

#SimulationPath
SimPath = "/Users/chiche/Desktop/DeepCrSearch/Simulations/DeepCrLib"
DataPath = "/Users/chiche/Desktop/DeepCrSearch/Analysis/SimData"
LibPath = "/Users/chiche/Desktop/REGevents"

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

    Ncores = 5
    IdAnt = np.random.randint(0, len(AntPow), Ncores)
    AntPow =  AntPow[IdAnt,:]
    Traces_C_deep = {k: Traces_C_pow[k] for k in IdAnt}
    Traces_G_deep = {k: Traces_G_pow[k] for k in IdAnt}

    for j in range(Ncores):
        
        SaveEvent(AntPow, Traces_C_deep, Traces_G_deep, Selfile, LibPath, j)
    
    sys.exit()
   
        
         
