#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 18:53:42 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt

# radio extent at 10^16.5 eV
re0316 = np.array([320, 352, 369, 393, 418, 410, 474, 546, 540])
# radio extent at 10^17 eV
re01 = np.array([578, 554, 590, 628, 670, 778, 834, 884])

# Corresponding zenith angles
theta0316 = np.array([0, 10, 20, 28, 34, 39 ,43 ,47, 50])
theta01 = np.array([10, 20, 28, 34, 39 ,43 ,47, 50])

# plot
plt.plot(theta0316, re0316)
# plot with cos(theta correction)
plt.plot(theta0316, re0316[0]/np.cos(theta0316*np.pi/180))
plt.show()

# plots with correction in (1+k*theta)/cos(theta)
plt.plot(theta0316, re0316)
plt.plot(theta0316, re0316[0]*(1+0.2062*theta0316*np.pi/180.0)/np.cos(theta0316*np.pi/180))
plt.show()

# correction in (1+k*theta)/cos(theta) at 10^17
plt.plot(theta01, re01)
plt.plot(theta01, re01[0]*(1+0.2062*theta01*np.pi/180.0)/np.cos(theta01*np.pi/180)/1.18)
plt.show()


