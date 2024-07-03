#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 17:06:13 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt

Altitude = 3216
Depth = np.array([100,80,60, 40, 0])
step = 10
limit = 150


def RectGrid(Altitude, Depth, step, limit):
    
    xp = np.arange(-limit, (limit + step) ,step)
    yp = np.arange(-limit, (limit + step) ,step)

    X, Y = np.meshgrid(xp, yp)
    xrect = X.ravel()
    yrect = Y.ravel()
    Nrect = len(xrect)
    
    xall, yall, zall = np.array([]), np.array([]), np.array([])
    
    for i in range(len(Depth)):
        
        xall = np.concatenate([xall, xrect])
        yall = np.concatenate([yall, yrect])
        z = Altitude - Depth[i]
        zall = np.concatenate([zall, np.ones(len(xrect))*z])

    return xall, yall, zall, Nrect

xall, yall, zall, Nrect =  RectGrid(Altitude, Depth, step, limit)



Plot = True

if(Plot):
    ax = plt.figure().add_subplot(projection='3d')
    ax.scatter(xall, yall, zs=zall, label='curve in (x, y)')    
    #plt.scatter(xall, yall)
    plt.savefig("/Users/chiche/Desktop/RectGridFull.pdf", bbox_inches = "tight")
    plt.show()

