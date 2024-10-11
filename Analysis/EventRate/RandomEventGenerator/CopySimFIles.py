#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Oct 10 10:52:23 2024

@author: chiche
"""


import numpy as np
import os 
import glob
import shutil

DataPath =  "/Users/chiche/Desktop/SimFiles/"
DataDir = "/Users/chiche/Desktop/DataDir/Atmosphere.dat"
Esim = glob.glob(DataPath + "*")

for i in range(len(Esim)):
    
    thetasim = glob.glob(Esim[i] + "/*") 
    for j in range(len(thetasim)):
        targetdir = thetasim[j] + "/"
        shutil.copy(DataDir, targetdir)
    
 