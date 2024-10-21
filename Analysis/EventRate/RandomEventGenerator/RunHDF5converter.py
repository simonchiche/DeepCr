#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Oct 12 15:02:08 2024

@author: chiche
"""

import os
import subprocess

def RunHDF5converter(HDF5path, REASpath):
    
    filemame = REASpath.split("/")[-1] + ".hdf5"
    
    cmd = "python " +  HDF5path + "/coreas_to_hdf5_airice.py --add_faerie_simulation " \
    + REASpath + "/SIM000001.reas -of " + filemame
    
    os.chdir(HDF5path)
    
    p=subprocess.Popen(cmd, cwd=os.getcwd(), shell=True)
    stdout, stderr = p.communicate()
    
    return 
