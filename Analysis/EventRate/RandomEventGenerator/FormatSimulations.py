#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 14:38:32 2024

@author: chiche
"""

import sys
import numpy as np
import matplotlib.pyplot as plt
from GetRandomDirection import GetRandomTheta, GetRandomPhi
from GenRandomE import GetRandomE
from CreateRandomEvent import GetDeepAntennaLayerEvents, SaveEvent
from SelectClosestEvent import LoadClosestEvent
from RunHDF5converter import RunHDF5converter
from GetFaerieWaveforms import GetFaerieVoltage
from icecream import ic
import glob

ic("Initalizing the paths")

# Simulation Path
SimPath = "/Users/chiche/Desktop/DeepCrSearch/Simulations/DeepCrLib"
DataPath = "/Users/chiche/Desktop/DeepCrSearch/Analysis/SimData/DeepCrLib"
DataFilesPath = "/Users/chiche/Desktop/REGevents/SimFiles/"
LibPath = "/Users/chiche/Desktop/REGevents/Events"
HDF5_LibPath = "/Users/chiche/Desktop/REGevents/HDF5files"
PythonPath = "/Users/chiche/Desktop/REGevents/coreas_to_hdf5_airice.py"

# ground level in meters
glevel =3216

# Number of events 
Nev = 50

# zenith limits
min_th = 0
max_th = 50
theta_rand = GetRandomTheta(min_th, max_th, Nev)

# azimuth limits 
min_phi = 0
max_phi = 360  
phi_rand = GetRandomPhi(min_phi, max_phi, Nev)

plt.hist(phi_rand)

# Energy limits (log)
pmin = 16.5
pmax = 17.5
#Eall = np.array([1e16, 3.16e16, 1e17, 3.16e17])
Eall = np.array([1e16])
PhiAll = np.array([0])
#ZenithAll =  np.array([0,10, 20, 28, 34, 39, 43, 47, 50])
ZenithAll =  np.array([34, 39, 43])
#GetRandomE(pmin, pmax, Nev)
print(Eall)
detectordescription = "/Users/chiche/Desktop/DeepCrSearch/NuRadioMC/NuRadioReco/detector/RNO_G/RNO_single_channel.json"

TriggerAllSims = dict()
AntposAll = dict()
EventPathAll = []
HDF5pathAll = []
#ThetaAll = dict() theta_rand
#EnergyAll = dict() energy_rand

GenerateEvent = True
if(GenerateEvent):
        for i in range(0,len(Eall)):
                for j in range(len(ZenithAll)):
                    print("\n")
                    ic('Event', i)
                    ic(Eall[i]/1e18, ZenithAll[j], PhiAll[0])

                    # Select the closest event and extract the traces for antennas at a depth of 100 meters
                    ic('Select closest event and get deep antennas traces')
                    ###AntPow, Traces_C_pow, Traces_G_pow, xrand, yrand, Selfile =  GetDeepAntennaLayerEvents(0.1, 10, 0, glevel, SimPath, DataPath)
                    #AntPow, Traces_C_pow, Traces_G_pow, xrand, yrand, Selfile =  GetDeepAntennaLayerEvents(Eall[i], ZenithAll[j], PhiAll[0], glevel, SimPath, DataPath)
                    Traces_C, Traces_G, AntPos, Selfile = \
                        LoadClosestEvent(Eall[i], ZenithAll[j], PhiAll[0], SimPath, DataPath)
                    print(Selfile)
                    
                    #break
                    # Save the event in the DataFilesPath in the format expected by the HDF5 converter in the EventPath repository
                    ic("saving event")
                    EventPath = LibPath + "/" + Selfile + "_" + str(0)
                    SaveEvent(AntPos, Traces_C, Traces_G, Selfile, DataFilesPath, EventPath)
                    ic("event saved")
                    #sys.exit()
                    EventPathAll.append(EventPath)

                    HDF5_FilePath = RunHDF5converter(PythonPath, EventPath, HDF5_LibPath)
                    ic("HDF5 conversion successful")

                    HDF5pathAll.append(HDF5_FilePath)
                    #print("savedir:", EventPath)
                    #print("HDF5path:", HDF5_LibPath)
                    #print("HDF5_FilePath:", HDF5_FilePath)
                    #sys.exit()