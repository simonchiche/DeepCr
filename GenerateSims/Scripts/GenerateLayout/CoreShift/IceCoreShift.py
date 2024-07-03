#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 20:11:21 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt
import sys 
sys.path.append('../Layouts/CrossGrid')
from CrossGrid import CrossGrid
from IceCorePosition import core

Altitude = 3216
Depth = np.array([100,80,60, 40])
step = 10
limit = 100
theta = 0
model = 1 # Greenland
phi = 0 #azimuth

xall, yall, zall, Ncross = CrossGrid(Altitude, Depth, step, limit)

def GetCoreShift(xall, yall, zall, Ncross, Depth, theta, phi, model):

    for i in range(len(Depth)):
        rcore = core(Depth[i], theta, model)[0][-1]
    
        xcore = -rcore*np.cos(phi*np.pi/180)
        ycore = -rcore*np.sin(phi*np.pi/180)
        print(xcore, ycore)
        xall[(i)*Ncross:(i+1)*Ncross] = xall[(i)*Ncross:(i+1)*Ncross] + xcore
        yall[(i)*Ncross:(i+1)*Ncross] = yall[(i)*Ncross:(i+1)*Ncross] + ycore
        
    return xall, yall, zall, Ncross

# =============================================================================
#                           Test Plot
# =============================================================================

xall, yall, zall, Ncross = \
GetCoreShift(xall, yall, zall, Ncross, Depth, theta, phi, model)

Plot = True
if(Plot):
    k = 3
    plt.scatter(xall[k*Ncross:(k+1)*Ncross], yall[k*Ncross:(k+1)*Ncross])
    plt.xlim(-10,10)

