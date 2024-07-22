#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 10:38:03 2024

@author: chiche
"""



import numpy as np


def GetDepthScaling(Depth):

    p0 = 31
    p1 = 2.813e-5
    p2 = 3.07

    r0 = p0

    r = p0 + p1*Depth**p2

    k = r/r0

    return k


def Getradius(x, A, B,C):
    
    # Function to fit the dependence of the radio footprint with depth
    # results to be indicated here in comments following "theta__: A = , B = , C = "
    #theta0: A = 31, B = 2.813e-5, C = 3.07
    #theta10: A = 319, B = 4.1152e-4, C = 3.032
    #theta20: A = 336, B = 1.4380e-3, C = 2.77
    #theta28: A = 376, B = 1.9917e-4, C = 3.23
    #theta34: A = 407, B = 5.4562e-2, C = 1.932
    #theta39: A = 403, B = 1.7811e-1, C = 1.703
    #theta43: A = 461, B = 2.5802e-1, C = 1.664
    #theta47: A = 534, B = 9.0471e-3, C = 2.446
    #theta50: A = 538, B = 6.1912e-4, C = 2.0719
    #p0_0316 = np.array([31, 319, 336, 376, 407, 403, 461, 534, 538])
    #p1_0316 = np.array([2.813e-5, 4.1152e-4, 1.4380e-3, 1.9917e-4, 5.4562e-2, 1.7811e-1, 2.5802e-1, 9.0471e-3, 6.1912e-2])
    #p2_0316 = np.array([3.07, 3.032, 2.77, 3.23, 1.932, 1.703, 1.664, 2.446, 2.0719])   
    # x: depth
    
    r = (A + B*x**C)
    
    return r


def GetEnergyScaling(E):

    E0 = 1e16
    k = (E/E0)**0.327

    return k