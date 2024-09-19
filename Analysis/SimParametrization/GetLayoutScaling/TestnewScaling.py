#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug 12 17:26:06 2024

@author: chiche
"""

import numpy as np
import matplotlib.pyplot as plt

from Modules.SimParam.GetLayoutScaling import GetDepthScaling, GetEnergyScaling, GetZenithScaling, ScalingFactor

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


def ScalingFactorOld(E, theta, Depth):
    
    k= GetDepthScalingOld(Depth)*GetEnergyScalingOld(E)*GetZenithScalingOld(theta)
    
    return k


ScalingFactor(1e16, 0, 0)
ScalingFactorOld(1e16, 0, 0)


theta_all = np.linspace(0,50, 20)
energy_all = np.array([1e16, 10**16.5, 10**17, 10**17.5, 10**18])
depth_all = np.array([0, 40, 60, 80, 100])


ktheta = np.zeros(len(theta_all))
ke  = np.zeros(len(energy_all))
kdepth =  np.zeros(len(depth_all))

kthetaOld = np.zeros(len(theta_all))
keOld = np.zeros(len(energy_all))
kdepthOld =  np.zeros(len(depth_all))

for i in range(len(ktheta)):
   
    ktheta[i] = ScalingFactor(1e16, theta_all[i], 0)
    kthetaOld[i] = ScalingFactorOld(1e16, theta_all[i], 0)

for i in range(len(ke)):
    ke[i] = ScalingFactor(energy_all[i], 0, 0)
    keOld[i] = ScalingFactorOld(energy_all[i], 0, 0)
    
for i in range(len(kdepth)):
    kdepth[i] = ScalingFactor(1e16, 0, depth_all[i])
    kdepthOld[i] = ScalingFactorOld(1e16, 0, depth_all[i])
      
    
plt.scatter(theta_all, ktheta)  
plt.scatter(theta_all, kthetaOld) 
plt.show()  

plt.scatter(energy_all, ke)  
plt.scatter(energy_all, keOld) 
plt.show()  

plt.scatter(depth_all, kdepth)  
plt.scatter(depth_all, kdepthOld) 
plt.show()  
 

    