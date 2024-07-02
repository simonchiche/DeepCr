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


def GetEnergyScaling(E):

    E0 = 1e16
    k = np.sqrt(E/E0)

    return k