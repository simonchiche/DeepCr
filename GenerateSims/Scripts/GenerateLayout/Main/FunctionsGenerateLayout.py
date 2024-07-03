#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 16:35:33 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt
from IceParameters import *
from GenerateStarshapes import generate_arm
import subprocess
import math
import os
import sys
import glob


def generate_input_file(x_pos, y_pos, z_pos, NantMax, SimDir):
    
    # Total number of antennas
    Nant = len(x_pos)
    
    # Number of group of antennas
    Ngroup = math.ceil(Nant/NantMax)
    
    Simulations = glob.glob(SimDir + "/*")
    
    # Loop over all the simulations of a given directory
    for j in range(len(Simulations)):
        filename = Simulations[j].split("/")[-1] #os.getcwd().split("/")[-1]
        
        # Loop over all the antenna groups a given antenna list
        for i in range(Ngroup):
        
            # Limits of the antennas number that will be simulated for this step
            start = i*NantMax
            stop = (i+1)*NantMax
            if(stop>Nant): stop = Nant
            
            # We create a new directory for eeach group and go into it
            path = Simulations[j] + "/" + filename  + "_"\
                + str(start) + "_" +  str(stop)
            cmd = "mkdir -p " + path
            #print(filename)
            #print(cmd)
            #sys.exit()
            p =subprocess.Popen(cmd, cwd=os.getcwd(), shell=True)
            stdout, stderr = p.communicate()
            
            os.chdir(path)
            
            file = open('SIM.list', 'w')
            
    
            # Loop over all the antennas of a given antenna group
            for k in range(start, stop, 1):
                file.write("AntennaPosition = %.1f %.1f %.1f ch%.d \n" \
                           %(x_pos[k], y_pos[k], z_pos[k], k))
                
            file.close()