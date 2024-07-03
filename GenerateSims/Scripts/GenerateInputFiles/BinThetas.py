#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 29 21:57:44 2024

@author: chiche
"""


import numpy as np 


min_th = np.log(np.cos(0))
max_th= np.log(np.cos(50*np.pi/180))

print(np.arccos(np.exp(np.linspace(min_th, max_th, 8)))*180/np.pi)


min_th = np.cos(0)
max_th= np.cos(50*np.pi/180)


print(np.arccos(np.linspace(min_th, max_th, 8))*180/np.pi)