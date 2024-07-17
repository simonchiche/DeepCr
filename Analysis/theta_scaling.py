#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Jul 16 18:53:42 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt


re0316 = np.array([320, 352, 369, 393, 418, 410, 474, 546, 540])

theta0316 = np.array([0, 10, 20, 28, 34, 39 ,43 ,47, 50])


plt.plot(theta0316, re0316)
plt.plot(theta0316, re0316[0]/np.cos(theta0316*np.pi/180))
plt.show()


plt.plot(theta0316, re0316)
plt.plot(theta0316, re0316[0]*(1+0.2062*theta0316*np.pi/180.0)/np.cos(theta0316*np.pi/180))
plt.show()


re01 = np.array([578, 554, 590, 628, 670, 778, 834, 884])
theta01 = np.array([10, 20, 28, 34, 39 ,43 ,47, 50])


plt.plot(theta01, re01)
plt.plot(theta01, re01[0]*(1+0.2062*theta01*np.pi/180.0)/np.cos(theta01*np.pi/180)/1.18)
plt.show()


def ScalingFactor(power, theta, k):
    
    k= (1+k*theta*np.pi/180.0)*(1e17/3.16e16)**power/(1e17/3.16e16)**0.5
    
    return k


ScalingFactor(0.36, 39, 0.2062)
