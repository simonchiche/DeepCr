#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Oct  4 00:42:23 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt

Bx = 7.7
Bz = -54.111

ub = np.array([Bx, 0, Bz])
ub_false = np.array([Bx, 0, -Bz])

theta = (180 - np.array([0,10,20,28,34,39, 43,47,50]))*np.pi/180.0

phi =0
sin_alpha = np.zeros(len(theta))
sin_alpha_false = np.zeros(len(theta))

for i in range(len(theta)):

    uv = np.array([np.sin(theta[i]) * np.cos(phi), np.sin(theta[i]) * np.sin(phi), np.cos(theta[i])])
    
    alpha = np.arccos(np.dot(uv, ub) / (np.linalg.norm(uv) * np.linalg.norm(ub)))
    alpha_false = np.arccos(np.dot(uv, ub_false) / (np.linalg.norm(uv) * np.linalg.norm(ub_false)))
    
    sin_alpha[i] = np.sin(alpha)
    sin_alpha_false[i] = np.sin(alpha_false)
    
    
plt.scatter(180-theta*180/np.pi, sin_alpha*180/np.pi)
plt.scatter(180-theta*180/np.pi, sin_alpha_false*180/np.pi)
plt.show()

plt.scatter(180-theta*180/np.pi, sin_alpha/sin_alpha_false)
plt.show()


# Projection du produit vectoriel
# simulation à 16.5 eV
# google doc + wiki

