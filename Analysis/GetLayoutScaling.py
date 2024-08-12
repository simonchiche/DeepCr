#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 20 10:38:03 2024

@author: chiche
"""


import numpy as np


def GetDepthScaling(Depth):

    p0 = 534
    p1 = 0.0090471
    p2 = 2.446

    r0 = p0

    r = p0 + p1*Depth**p2

    k = r/r0

    return k



def GetEnergyScaling(E):

    E0 = 1e16
    k = (E/E0)**0.327

    return k

def GetZenithScaling(theta):

    k= (1 + 0.308*theta*np.pi/180.0)/np.cos(theta*np.pi/180.0)
    return k


# =============================================================================
#                               Tests
# =============================================================================

def ScalingFactor(E, theta, Depth):
    
    k= GetDepthScaling(Depth)*GetEnergyScaling(E)*GetZenithScaling(theta)
    
    return k

