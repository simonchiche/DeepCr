#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 15:02:08 2024

@author: chiche
"""

import os
import subprocess
from icecream import ic

def RunHDF5converter(PythonPath, REASpath, HDF5path):

    """"
    Create an hdf5 file, from the path to a folder containing a reas file and formatted as expected by the converter
    HDF5path: path where the file is stored
    REASpath: path to the REAS file
    """
    os.chdir(REASpath)
    ic('HDF5 conversion')
    #os.chdir(REASpath)
    filename = REASpath.split("/")[-1] + ".hdf5"
    
    cmd = "python " +  PythonPath  + " --add_faerie_simulation " \
    + REASpath + "/SIM000001.reas -of " + HDF5path + "/" + filename
    
    ic(cmd)
    print(cmd)
    
    p=subprocess.Popen(cmd, cwd=os.getcwd(), shell=True)
    stdout, stderr = p.communicate()

    HDF5_FilePath = HDF5path + "/" + filename
    ('HDf5 conversion successful')
    return HDF5_FilePath

####HDF5path = "/Users/chiche/Desktop/Test"
#####REASpath = HDF5path + "/Rectangle_Proton_0.001_0_0_1"
'''
PythonPath = "/Users/chiche/Desktop/REGevents/coreas_to_hdf5_airice.py"
REASpath = "/Users/chiche/Desktop/REGevents/Events/Rectangle_Proton_0.0316_34_0_1_2/" #"/Users/chiche/Desktop/REGevents/Events/Rectangle_Proton_0.1_20_0_1_0" 
HDF5path =  "/Users/chiche/Desktop/REGevents/HDF5files" #REASpath  #"/Users/chiche/Desktop/HDF5_converte_test/Rectangle_Proton_0.001_0_0_1"


RunHDF5converter(PythonPath, REASpath, HDF5path)
#ic("saved")
'''