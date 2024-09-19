#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 11:33:42 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt


def InAirRayTracing(Xmax, theta, phi):
    
    #------------------------------#
    # Xmax: 3D numpy array (x_Xmax, y_Xmax, z_Xmax)
    
    #------------------------------#
    
    # Get shower direction
    uv = ...
    
   # aperture angle of the cone of radio emission
   ConeAngle = 5 
   # resolution angle of the ray-tracing procedure
   ResAngle = 0.1 
   # Number of rays
   Nrays = int(ConeAngle/ResAngle) +1
   
   LaunchAngles =  np.linspace(0, 5, Nrays)
   
   for i in range(len(LaunchAngles)):
       
       # Propagate each ray from Xmax to the ground
       
    