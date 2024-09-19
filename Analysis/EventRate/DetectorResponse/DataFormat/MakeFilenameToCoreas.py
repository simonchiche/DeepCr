#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Aug 20 11:40:30 2024

@author: chiche
"""


import numpy as np
import subprocess
import glob
import os
import sys
path = "/Users/chiche/Desktop/Rectangle_Proton_0.1_20_0_1/SIM000001_geant"

GeantFiles = glob.glob(path + "/*")

for i in range(len(GeantFiles)):
    
    filename_root = GeantFiles[i].split("/")[-1].split(".")[0]
    antnum = filename_root.split("a")[1]
    
    cmd = "mv " + GeantFiles[i] + " " + path + "/raw_ch" + antnum + ".dat"
    p = subprocess.Popen(cmd, cwd=os.getcwd(), shell=True)
    stdout, stderr = p.communicate()
    #if(i<10): print(cmd)