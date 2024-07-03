#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Mar 14 11:37:12 2024

@author: chiche
"""


import numpy as np
import subprocess
import matplotlib.pyplot as plt
import math
from IceParameters import *
#print(IceParameters.step)

#IceLevel = 2835*1e2
#AntennaDepth = -100*1e2


# =============================================================================
#                        Layout Building
# =============================================================================

step = 20*1e2
AntArm = 20
angles = np.array([0, 45, 90, 135, 180, 225, 270, 315])

def generate_arm(AntArm, step, angles, IceLevel, AntennaDepth):
    
    x_pos, y_pos, z_pos = [], [], []
    
    for i in range(AntArm):
        for j in range(len(angles)):
        
            x = (i+1)*step*np.cos(angles[j]*np.pi/180.0)
            y = (i+1)*step*np.sin(angles[j]*np.pi/180.0)
    
            x_pos.append(x)
            y_pos.append(y)

    z_pos = (IceLevel + AntennaDepth)*np.ones(len(x_pos))
    
    return np.array(x_pos),np.array(y_pos), np.array(z_pos)

# We generate the antennas positions
x_pos, y_pos, z_pos = \
generate_arm(AntArm, step, angles, IceLevel, AntennaDepth)

# We plot the result as a check
plt.scatter(x_pos/1e2, y_pos/1e2)
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.show()

