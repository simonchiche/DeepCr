#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 02:21:35 2024

@author: chiche
"""
#region Modules
import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt
import glob
import sys
import pickle
sys.path.append("/Users/chiche/Desktop/DeepCrSearch/Analysis/")
from Modules.SimParam.GetLayoutScaling import GetDepthScaling, GetEnergyScaling
from scipy.interpolate import interp1d
import scipy
from Modules.Fluence.FunctionsGetFluence import  Norm, LoadTraces, GetPeakTraces, Traces_cgs_to_si, GetDepths, CorrectScaling, CombineTraces, CorrectLength, GetIntTraces, GetIntTracesSum, GetRadioExtent
from Modules.Fluence.FunctionsPlotFluence import EfieldMap, PlotLDF, PlotTraces, plot_polarisation, PlotMaxTraces, PlotAllTraces, PlotLayer, PlotGivenTrace, PlotAllChannels
from CleanCoreasTraces import CleanCoreasTraces
from Modules.SimParam.PlotRadioSimExtent import PlotFillingFactor, PlotRadioSimExtent
from scipy.interpolate import griddata
from datetime import datetime
from scipy.optimize import curve_fit
from scipy.signal import butter, filtfilt
from Modules.Fluence.FunctionsRadiationEnergy import GetFluence, GetRadiationEnergy
#endregion

#region Plot settings
font = {'family' : 'DejaVu Sans',
        'weight' : 'normal',
        'size'   : 14}
plt.rc('font', **font)
#endregion

ZenithAll = np.array([0,10,20,28,34,39,43,47,50])
EnergyAll = np.array([0.0316, 0.1, 0.316])
DepthAll = np.array([0,40, 60, 80, 100])
glevel = 3216

#region Path definition
PowerDataPath = "/Users/chiche/Desktop/DeepCrSearch/Analysis/GetFluence/Data/Power/"
SimDir =     "DeepCrLib" #"InterpSim"
RootName = "Rectangle_Proton_" 
TailName =  "_0_1"
#endregion

header = "# Zenith\tEnergy\tEfield_x\tEfield_y\tEfield_z\tEfield_tot\tfx\tfy\tfz\tftot"

with open(PowerDataPath + "AirPowerMap.txt", "a") as f1, open(PowerDataPath + "IcePowerMap.txt", "a") as f2:
    #f1.write(header + "\n")
    #f2.write(header + "\n")
    for i in range(len(EnergyAll)):
        for j in range(len(ZenithAll)):
            for k in range(len(DepthAll)):
                print(i,j,k)
                SimName = RootName + str(EnergyAll[i]) + "_" + str(ZenithAll[j]) + TailName
                zenith = ZenithAll[j]
                
                
                # region Data loading and formatting
                SimDataPath = "/Users/chiche/Desktop/DeepCrSearch/Analysis/SimData/" + SimDir + "/" + SimName 

                Nant  = np.load(SimDataPath + "/Nant.npy")
                with open(SimDataPath + '/Traces_C.pkl', 'rb') as file:
                    Traces_C = pickle.load(file)
                with open(SimDataPath + '/Traces_G.pkl', 'rb') as file:
                    Traces_G = pickle.load(file)
                Pos  = np.load(SimDataPath + "/Pos.npy", allow_pickle=True)

                with open(SimDataPath + '/Traces_tot.pkl', 'rb') as file:
                    Traces_tot = pickle.load(file)

                #cgs to si
                Traces_C = Traces_cgs_to_si(Traces_C)
                Traces_G = Traces_cgs_to_si(Traces_G)
                #endregion

                # =============================================================================
                #                                 Get integral
                # =============================================================================

                ExC_int, EyC_int, EzC_int, EtotC_int = GetIntTraces(Traces_C, Nant)

                ExG_int, EyG_int, EzG_int, EtotG_int = GetIntTraces(Traces_G, Nant)

                Ex_tot_int, Ey_tot_int, Ez_tot_int, Etot_int = GetIntTraces(Traces_tot, Nant)

                # =============================================================================
                #                                 Write Power
                # =============================================================================

                zice = glevel - DepthAll
                spacing = abs(Pos[Pos[:,2]==zice[k]][1,0]- Pos[Pos[:,2]==zice[k]][0,0])
                print(spacing)
                #In-air emission
                #Efield
                EfieldAir_x, EfieldAir_y, EfieldAir_z, EfieldAir_tot = \
                    np.sum(ExC_int[Pos[:,2] == zice[k]]), np.sum(EyC_int[Pos[:,2] == zice[k]]),\
                    np.sum(EzC_int[Pos[:,2] == zice[k]]), np.sum(EtotC_int[Pos[:,2] == zice[k]])
                
                #Fluence
                fAir_x, fAir_y, fAir_z, fAir_tot = \
                    np.sum(ExC_int[Pos[:,2] == zice[k]]**2)*spacing**2, np.sum(EyC_int[Pos[:,2] == zice[k]]**2)*spacing**2,\
                    np.sum(EzC_int[Pos[:,2] == zice[k]]**2)*spacing**2, np.sum(EtotC_int[Pos[:,2] == zice[k]]**2)*spacing**2
                
                f1.write(f"{zenith}\t{EnergyAll[i]}\t{zice[k]}\t{EfieldAir_x:.3e}\t{EfieldAir_y:.3e}\t{EfieldAir_z:.3e}\t{EfieldAir_tot:.3e}\t"
                        f"\t{fAir_x:.3e}\t{fAir_y:.3e}\t{fAir_z:.3e}\t{fAir_tot:.3e}\n") 

                #In-ice emission
                #Efield
                EfieldIce_x, EfieldIce_y, EfieldIce_z, EfieldIce_tot = \
                    np.sum(ExG_int[Pos[:,2] == zice[k]]), np.sum(EyG_int[Pos[:,2] == zice[k]]),\
                    np.sum(EzG_int[Pos[:,2] == zice[k]]), np.sum(EtotG_int[Pos[:,2] == zice[k]])

                #Fluence
                fIce_x, fIce_y, fIce_z, fIce_tot = \
                    np.sum(ExG_int[Pos[:,2] == zice[k]]**2)*spacing**2, np.sum(EyG_int[Pos[:,2] == zice[k]]**2)*spacing**2,\
                    np.sum(EzG_int[Pos[:,2] == zice[k]]**2)*spacing**2, np.sum(EtotG_int[Pos[:,2] == zice[k]]**2)*spacing**2
                f2.write(f"{zenith}\t{EnergyAll[i]}\t{zice[k]}\t{EfieldIce_x:.3e}\t{EfieldIce_y:.3e}\t{EfieldIce_z:.3e}\t{EfieldIce_tot:.3e}\t"
                        f"\t{fIce_x:.3e}\t{fIce_y:.3e}\t{fIce_z:.3e}\t{fIce_tot:.3e}\n")   


                

