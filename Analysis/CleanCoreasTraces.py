#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 15:54:59 2024

@author: chiche
"""

import numpy as np

from FunctionsGetFluence import Norm
import matplotlib.pyplot as plt


def CleanCoreasTraces(Traces):
    
    window = 2500 # time window of the signal in ns
    
    Nant = len(Traces)
    print(Nant)
    CleanTraces = dict()
    for i in range(Nant):
        
        trace = Norm(Traces[i][:,1], Traces[i][:,2], Traces[i][:,3])
        
        binning = Traces[i][1,0] - Traces[i][0,0]
        binning =  round(binning*int(1e9), 1)
        
        imax = np.argmax(trace)
        ilow = int(imax -(window/(2.0*binning)))
        ihigh = int(imax + (window/(2.0*binning)))
        
        if(ilow<0):
            ilow = 0
            ihigh = int(window/binning)
        
        cleanT = Traces[i][ilow:ihigh, 0]
        cleanx = Traces[i][ilow:ihigh, 1]
        cleany = Traces[i][ilow:ihigh, 2]
        cleanz = Traces[i][ilow:ihigh, 3]
        
        List = [cleanT, cleanx, cleany, cleanz]
        CleanTraces[i] = np.array(List).T
        
    return CleanTraces
        
        

