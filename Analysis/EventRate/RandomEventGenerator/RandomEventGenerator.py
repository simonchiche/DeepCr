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

# Number of events 
Nev = 5000
# zenith limits
min_th = 0
max_th = 50
# azimuth limits 
min_phi = 0
max_phi = 360  

theta_rand = GetRandomTheta(min_th, max_th, Nev)
phi_rand = GetRandomPhi(min_phi, max_phi, Nev)