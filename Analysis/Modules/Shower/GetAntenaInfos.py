#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 13:29:49 2024

@author: chiche
"""

import numpy as np

def GetDepths(AntPos, glevel):
    
    zpos = AntPos[:,2]
    Depths, indices = np.unique(zpos, return_index=True)
    Depths = Depths[np.argsort(indices)]
    Depths = abs(Depths-glevel)
    Nlay = len(Depths)
    NantLay = int(len(AntPos)/Nlay)
    
    return Nlay, NantLay, Depths

def GetAntennaLayer(AntPos, glevel, Depth):
    
    targetz = glevel - Depth
    AntennaLayer =  AntPos[AntPos[:,2] == targetz]
    
    return AntennaLayer

def GetSurfaceAntennas(AntPos, glevel):
    
    SurfaceAntennas = GetAntennaLayer(AntPos, glevel, 0)
    
    return SurfaceAntennas

