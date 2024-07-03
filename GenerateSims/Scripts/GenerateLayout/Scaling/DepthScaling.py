#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 13 10:59:41 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit



depth = np.array([100,80,60,40])
radius = np.array([70,50,40,33])


plt.scatter(depth, radius)


def Getradius2(x, A, B,C):
    
    # x : depth
    
    r = A + B*x**C
    
    return r

def Getradius(x, C):
    
    # x : depth
    
    r60 = 40
    r = r60*(C+x)/(C+60)
    
    return r



popt, pcov = curve_fit(Getradius2, depth, radius)


k60 =  Getradius2(0, popt[0], popt[1], popt[2])\
/Getradius2(60, popt[0], popt[1], popt[2])
print(k60)

def GetDepthScaling(Depth):
    
    p0 = 31
    p1 = 2.813e-5
    p2 = 3.07
    
    r0 = p0
    
    r = p0 + p1*Depth**p2
    
    k = r/r0
    
    return k
    
print(GetDepthScaling(40))
