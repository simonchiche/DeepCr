#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 23 14:05:18 2024

@author: chiche
"""

import numpy as np

theta = 10*180/np.pi
phi =np.array([0,90,180])*180/np.pi
i = 75*180/np.pi

for k in range(len(phi)):

    uv = np.array([np.sin(theta)*np.cos(phi[k]), np.sin(theta)*np.sin(phi[k]), np.cos(theta)])
    
    Bx =  np.cos(i)
    Bz = - np.sin(i)
    
    uB = np.array([Bx, 0, Bz])
    
    
    uv_x_uB = np.cross(uv, uB) # unit vector along the vxb direction
    uv_x_uB /= np.linalg.norm(uv_x_uB) # normalisation
    
    uv_x_uvxB  = np.cross(uv, uv_x_uB) # unit vector along the vxvxb direction
    uv_x_uvxB /= np.linalg.norm(uv_x_uB) # normalisation
    
    
    P = np.transpose(np.array([uv, uv_x_uB, uv_x_uvxB])) 
    # matrix to go from the shower reference frame to the geographic reference frame
        
    P_inv = np.linalg.inv(P) 
    print(P_inv)
    print("\n")
    

