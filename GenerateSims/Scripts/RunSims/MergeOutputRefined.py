#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr 18 17:50:52 2024

@author: chiche
"""



import glob
import subprocess
import os

HostDir = "CoreShiftTest"
OutputPath = "/pnfs/iihe/rno-g/user/schiche/DeepCr/AirIceSims/OutputFiles/FinalOutput/"
OutputPath = OutputPath + HostDir 

SimPath = glob.glob(OutputPath + "/*")
Nsim = len(SimPath)

for k in range(Nsim):

    GeantFiles = glob.glob(SimPath[k] + "/Geant/*.dat")
    CoreasFiles =  glob.glob(SimPath[k] + "/Coreas/*.dat")
#GeantFiles = glob.glob("/pnfs/iihe/rno-g/user/schiche/DeepCr/AirIceSims/OutputFiles/FinalOutput/CoreShiftTest/Cross_Proton_0.01_0_0_1/Geant/*.dat")
#CoreasFiles = glob.glob("/pnfs/iihe/rno-g/user/schiche/DeepCr/AirIceSims/OutputFiles/FinalOutput/CoreShiftTest/Cross_Proton_0.01_0_0_1/Coreas/*.dat")

    # Create a list to store file paths
    Geantfile_paths = []
    Coreasfile_paths = []
    
    # Add file paths to the list
    for i in range(len(GeantFiles)):
        Geantfilename =  "antenna" + str(i) + ".dat"
        Geantfile_paths.append(Geantfilename)
    
        Coreasfilename =  "raw_ch" + str(i) + ".dat"
        Coreasfile_paths.append(Coreasfilename)

print(Coreasfile_paths)

    # Concatenate files using cat command
    filename = SimPath[k].sp;it("/"[-1])
    os.chdir(SimPath[k] + "/Geant/")
    subprocess.call(["cat"] + Geantfile_paths, stdout=open(filename + "_Geant.dat", "w"))
    
    
    # Concatenate files using cat command
    os.chdir(SimPath[k] + "/Coreas/")
    subprocess.call(["cat"] + Coreasfile_paths, stdout=open(filename + "_Coreas.dat", "w"))
    

