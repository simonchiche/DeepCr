#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 21:52:41 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt

def Getradius(x, A, B,C):
    
    # Function to fit the dependence of the radio footprint with depth    
    r = (A + B*x**C)
    
    return r

Depths = np.array([0,40,60,80,100]) + 3116


### Parametrization of the radio extent vs depth for different
# zenith angles and energies (to be saved in specific files) 
# see conclusions below

#fit results
p0_0316 = np.array([31, 319, 336, 376, 407, 403, 461, 534, 538])
p1_0316 = np.array([2.813e-5, 4.1152e-4, 1.4380e-3, 1.9917e-4, 5.4562e-2,\
                    1.7811e-1, 2.5802e-1, 9.0471e-3, 6.1912e-2])
p2_0316 = np.array([3.07, 3.032, 2.77, 3.23, 1.932, 1.703, 1.664, 2.446, 2.0719])   

p0_01 = np.array([555, 553, 582, 645, 667, 789, 823, 922])
p1_01 = np.array([3.5956e-2, 8.0731e-2, 1.3556e-1, 9.4164e-2, 3.5645e-1, \
                  1.4066e-2, 3.7933e-1, 4.4415e-2])
p2_01 = np.array([2.04, 1.89, 1.842, 1.911, 1.691, 2.45, 1.7652, 2.27])   


r0 = Getradius(abs(np.array(Depths)-3216), p0_0316[0], p1_0316[0], p2_0316[0])
r1 = Getradius(abs(np.array(Depths)-3216), p0_0316[1], p1_0316[1], p2_0316[1])
r2 = Getradius(abs(np.array(Depths)-3216), p0_0316[2], p1_0316[2], p2_0316[2])
r3 = Getradius(abs(np.array(Depths)-3216), p0_0316[3], p1_0316[3], p2_0316[3])
r4 = Getradius(abs(np.array(Depths)-3216), p0_0316[4], p1_0316[4], p2_0316[4])
r5 = Getradius(abs(np.array(Depths)-3216), p0_0316[5], p1_0316[5], p2_0316[5])
r6 = Getradius(abs(np.array(Depths)-3216), p0_0316[6], p1_0316[6], p2_0316[6])
r7 = Getradius(abs(np.array(Depths)-3216), p0_0316[7], p1_0316[7], p2_0316[7])
r8 = Getradius(abs(np.array(Depths)-3216), p0_0316[8], p1_0316[8], p2_0316[8])

r0_17 = Getradius(abs(np.array(Depths)-3216), p0_01[0], p1_01[0], p2_01[0])
r1_17 = Getradius(abs(np.array(Depths)-3216), p0_01[1], p1_01[1], p2_01[1])
r2_17 = Getradius(abs(np.array(Depths)-3216), p0_01[2], p1_01[2], p2_01[2])
r3_17 = Getradius(abs(np.array(Depths)-3216), p0_01[3], p1_01[3], p2_01[3])
r4_17 = Getradius(abs(np.array(Depths)-3216), p0_01[4], p1_01[4], p2_01[4])
r5_17 = Getradius(abs(np.array(Depths)-3216), p0_01[5], p1_01[5], p2_01[5])
r6_17 = Getradius(abs(np.array(Depths)-3216), p0_01[6], p1_01[6], p2_01[6])
r7_17 = Getradius(abs(np.array(Depths)-3216), p0_01[7], p1_01[7], p2_01[7])


### Plots
# At 10^16.5 eV
plt.plot(abs(np.array(Depths)-3216), r1/p0_0316[1])
plt.plot(abs(np.array(Depths)-3216), r2/p0_0316[2])
plt.plot(abs(np.array(Depths)-3216), r3/p0_0316[3])
plt.plot(abs(np.array(Depths)-3216), r0/p0_0316[0])

plt.plot(abs(np.array(Depths)-3216), r4/p0_0316[4])
plt.plot(abs(np.array(Depths)-3216), r5/p0_0316[5])
plt.plot(abs(np.array(Depths)-3216), r6/p0_0316[6])
plt.plot(abs(np.array(Depths)-3216), r7/p0_0316[7])
plt.plot(abs(np.array(Depths)-3216), r8/p0_0316[8])
plt.xlabel("Depth [m]")
plt.ylabel("Scaling factor")
plt.title("E = $10^{16.5}$ eV", fontsize =14)
#plt.savefig(OutputPath + "Scaling_vs_Depth_E%.2f.pdf" \
#             %(energy), bbox_inches = "tight")
plt.show()

# At 10^17 eV
plt.plot(abs(np.array(Depths)-3216), r0_17/p0_01[0])
plt.plot(abs(np.array(Depths)-3216), r1_17/p0_01[1])
plt.plot(abs(np.array(Depths)-3216), r2_17/p0_01[2])
plt.plot(abs(np.array(Depths)-3216), r3_17/p0_01[3])
plt.plot(abs(np.array(Depths)-3216), r4_17/p0_01[4])
plt.plot(abs(np.array(Depths)-3216), r5_17/p0_01[5])
plt.plot(abs(np.array(Depths)-3216), r6_17/p0_01[6])
plt.plot(abs(np.array(Depths)-3216), r7_17/p0_01[7])
plt.xlabel("Depth [m]")
plt.ylabel("Scaling factor")
plt.title("E = $10^{17}$ eV", fontsize =14)
#plt.savefig(OutputPath + "Scaling_vs_Depth_E%.2f.pdf" \
#             %(energy), bbox_inches = "tight")
plt.show()

# 10^16.5 eV and 10^17 eV mixed

plt.plot(abs(np.array(Depths)-3216), r1/p0_0316[1], linestyle =":")
plt.plot(abs(np.array(Depths)-3216), r2/p0_0316[2], linestyle =":")
plt.plot(abs(np.array(Depths)-3216), r3/p0_0316[3], linestyle =":")
plt.plot(abs(np.array(Depths)-3216), r0/p0_0316[0], linestyle =":")

plt.plot(abs(np.array(Depths)-3216), r4/p0_0316[4], linestyle =":")
plt.plot(abs(np.array(Depths)-3216), r5/p0_0316[5], linestyle =":")
plt.plot(abs(np.array(Depths)-3216), r6/p0_0316[6], linestyle =":")
plt.plot(abs(np.array(Depths)-3216), r7/p0_0316[7], linestyle =":")
plt.plot(abs(np.array(Depths)-3216), r8/p0_0316[8], linestyle =":")

plt.plot(abs(np.array(Depths)-3216), r0_17/p0_01[0])

plt.plot(abs(np.array(Depths)-3216), r1_17/p0_01[1])
plt.plot(abs(np.array(Depths)-3216), r2_17/p0_01[2])
plt.plot(abs(np.array(Depths)-3216), r3_17/p0_01[3])
plt.plot(abs(np.array(Depths)-3216), r4_17/p0_01[4])
plt.plot(abs(np.array(Depths)-3216), r5_17/p0_01[5])
plt.plot(abs(np.array(Depths)-3216), r6_17/p0_01[6])
plt.plot(abs(np.array(Depths)-3216), r7_17/p0_01[7])
plt.show()


# CCL: the seventh parametrization (in index) at 10^16.5 eV seems to reproduce well the dependency 
# of the radio extent with the depth independently of the simulation energy and depth
