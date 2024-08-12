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

def GetDepthScalingOld(Depth):

    p0 = 31
    p1 = 2.813e-5
    p2 = 3.07

    r0 = p0

    r = p0 + p1*Depth**p2

    k = r/r0

    return k

def GetEnergyScalingOld(E):

    E0 = 1e16
    k = (E/E0)**0.5

    return k

def GetZenithScalingOld(theta):

    k= 1/np.cos(theta*np.pi/180.0)
    return k

# =============================================================================
#                               Tests
# =============================================================================

def ScalingFactor(E, theta, Depth):
    
    k= GetDepthScaling(Depth)*GetEnergyScaling(E)*GetZenithScaling(theta)
    
    return k



ScalingFactor(1e16, 0, 0)


theta_all = np.linspace(0,50, 20)
energy_all = np.array([1e16, 10**16.5, 10**17, 10**17.5, 10**18])
depth_all = np.array([0, 40, 60, 80, 100])



