#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr  5 18:44:28 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt

Altitude = 3216
Depth = np.array([100,80,60, 40])
step = 10
limit = 100

def CrossGrid(Altitude, Depth, step, limit):

    minlim  = -limit
    
    xline = np.arange(minlim, (limit+ step), step)
    yline = np.zeros(len(xline))
    
    ycross = np.concatenate([yline, xline])
    xcross  = np.concatenate([xline, yline])
    Ncross = len(xcross)

     #= xcross, ycross
    xall, yall, zall = np.array([]), np.array([]), np.array([])

    for i in range(len(Depth)):
        
    
        xall = np.concatenate([xall, xcross])
        yall = np.concatenate([yall, ycross])
        z = Altitude - Depth[i]
        zall = np.concatenate([zall, np.ones(len(xcross))*z])
        
    return xall, yall, zall, Ncross

# =============================================================================
#                          Plot Test
# =============================================================================

xall, yall, zall, Ncross =  CrossGrid(Altitude, Depth, step, limit)

Plot = False

if(Plot):
    ax = plt.figure().add_subplot(projection='3d')
    ax.scatter(xall, yall, zs=zall, label='curve in (x, y)')    
    plt.scatter(xall, yall)
    plt.show()
