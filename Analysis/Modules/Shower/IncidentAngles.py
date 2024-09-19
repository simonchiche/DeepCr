#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Jul 22 15:54:02 2024

@author: chiche
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.append("/Users/chiche/Desktop/DeepCrSearch/Analysis/")
import Modules.Shower.GetXmaxCoordinate

def GetIncidentAngles(XmaxPos, AntPos):
    
    theta = np.zeros(len(AntPos))
    phi = np.zeros(len(AntPos))
    print(XmaxPos)
    for i in range(len(AntPos)):
        AntVec = AntPos[i]- XmaxPos
        
        Norm = np.sqrt(AntVec[0]**2 + AntVec[1]**2 + AntVec[2]**2) 
        uant = AntVec/Norm
        
        theta[i] = np.arccos(-uant[2])
        phi[i] = np.arctan2(uant[1], uant[0])
    
    return theta*180/np.pi, phi*180/np.pi

#theta_i, phi_i = GetIncidentAngles(XmaxPos, Pos)
