#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 14:38:32 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt
from GetRandomDirection import GetRandomTheta, GetRandomPhi
from GenRandomE import GetRandomE

#SimulationPath
SimPath = "/Users/chiche/Desktop/DeepCrSearch/Simulations/DeepCrLib"
DataPath = "/Users/chiche/Desktop/DeepCrSearch/Analysis/SimData"

# Number of events 
Nev = 5000

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


