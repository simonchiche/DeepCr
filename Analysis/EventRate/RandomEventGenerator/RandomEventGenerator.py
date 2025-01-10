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
from RunHDF5converter import RunHDF5converter
from GetFaerieWaveforms import GetFaerieVoltage
from icecream import ic
import glob
#ic.configureOutput(prefix=None)
ic("Initalizing the paths")

# Simulation Path
SimPath = "/Users/chiche/Desktop/DeepCrSearch/Simulations/DeepCrLib"
DataPath = "/Users/chiche/Desktop/DeepCrSearch/Analysis/SimData/DeepCrLib"
DataFilesPath = "/Users/chiche/Desktop/REGevents/SimFiles/"
LibPath = "/Users/chiche/Desktop/REGevents/Events"
HDF5_LibPath = "/Users/chiche/Desktop/REGevents/HDF5files"
PythonPath = "/Users/chiche/Desktop/REGevents/coreas_to_hdf5_airice.py"

ic("Initalizing the input parameters")

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
Erand = GetRandomE(pmin, pmax, Nev)

detectordescription = "/Users/chiche/Desktop/DeepCrSearch/NuRadioMC/NuRadioReco/detector/RNO_G/RNO_single_channel.json"

TriggerAllSims = dict()
AntposAll = dict()
EventPathAll = []
HDF5pathAll = []
#ThetaAll = dict() theta_rand
#EnergyAll = dict() energy_rand

GenerateEvent = True
if(GenerateEvent):
        for i in range(len(Erand)):
                print("\n")
                ic('Event', i)
                ic(Erand[i]/1e18, theta_rand[i], phi_rand[i])

                # Select the closest event and extract the traces for antennas at a depth of 100 meters
                ic('Select closest event and get deep antennas traces')
                ###AntPow, Traces_C_pow, Traces_G_pow, xrand, yrand, Selfile =  GetDeepAntennaLayerEvents(0.1, 10, 0, glevel, SimPath, DataPath)
                AntPow, Traces_C_pow, Traces_G_pow, xrand, yrand, Selfile =  GetDeepAntennaLayerEvents(Erand[i], theta_rand[i], phi_rand[i], glevel, SimPath, DataPath)
                break
                # Save the event in the DataFilesPath in the format expected by the HDF5 converter in the EventPath repository
                ic("saving event")
                EventPath = LibPath + "/" + Selfile + "_" + str(i)
                SaveEvent(AntPow, Traces_C_pow, Traces_G_pow, Selfile, DataFilesPath, EventPath)
                ic("event saved")
                EventPathAll.append(EventPath)

                HDF5_FilePath = RunHDF5converter(PythonPath, EventPath, HDF5_LibPath)
                ic("HDF5 conversion successful")

                HDF5pathAll.append(HDF5_FilePath)
                #print("savedir:", EventPath)
                #print("HDF5path:", HDF5_LibPath)
                #print("HDF5_FilePath:", HDF5_FilePath)
                #sys.exit()

GetTrigger = True
HDF5_FilePath = glob.glob("/Users/chiche/Desktop/REGevents/HDF5files/*")
ZenithAll = []
EnergyAll = []
PhiAll = []
if(GetTrigger):
        for i in range(len(Erand)):
                
                inputfilename = []
                inputfilename.append(HDF5_FilePath[i])

                VoltVpole, tVpole, Trigger, threshold,  freqChannel, SpectrumChannel, Antpos\
                = GetFaerieVoltage(inputfilename, detectordescription)
                ic('FAERIE reader succesful')
                TriggerAllSims[i] = Trigger
                AntposAll[i] = Antpos
                
                ZenithAll.append(int(HDF5_FilePath[i].split("/")[-1].split("_")[3]))
                EnergyAll.append(int(HDF5_FilePath[i].split("/")[-1].split("_")[2]))
                PhiAll.append(int(HDF5_FilePath[i].split("/")[-1].split("_")[4]))
                print(sum(Trigger), len(Trigger), sum(Trigger)/len(Trigger))
                print(len(VoltVpole))
                print(max(abs(VoltVpole[0])))  
                #sys.exit()
                '''
                #Threshold = 1e-3
                k =0
                ilist = np.array([203, 237, 264, 265])
                for i in range(len(tVpole)):
                #if(max(abs(VoltVpole[i])>Threshold) & (k<40)):
                #if(i in ilist):
                        print(i)
                        plt.plot(tVpole[i], VoltVpole[i])
                        plt.title("%.d" %i)
                        plt.xlabel('time / ns')
                        plt.ylabel('voltage / mV')
                        #plt.savefig("/Users/chiche/Desktop/DeepCrPlots_13_11_24/VoltageNoisy/VoltageTrace_%.d.pdf" %i, bbox_inches = "tight")
                        plt.show()
                        plt.plot(freqChannel[i], SpectrumChannel[i])
                        plt.xlabel('frequencies / GHz')
                        plt.title("%.d" %i)
                        #plt.savefig("/Users/chiche/Desktop/DeepCrPlots_13_11_24/VoltageNoisy/Spectrum_%.d.pdf" %i, bbox_inches = "tight")
                        plt.show()
                        k = k + 1

                ### Additionnal prints
                print(Trigger)
                print(len(VoltVpole))
                print(max(abs(VoltVpole[0])))

                plt.scatter(Antpos[:,0], Antpos[:,1], label = "All events")
                plt.scatter(Antpos[:,0][Trigger], Antpos[:,1][Trigger], color ="red", label = "Triggered events")
                plt.xlabel("x [m]")
                plt.ylabel("y [m]")
                plt.legend()
                plt.savefig("/Users/chiche/Desktop/Trigger.pdf", bbox_inches ="tight")
                plt.show()
         
                '''
'''            
ZenithSortCut = np.unique(ZenithAll)           
TriggRate = np.zeros(50)
MeanTriggRate = np.zeros(len(ZenithSortCut))
StdTriggRate = np.zeros(len(ZenithSortCut))

for i in range(50):
    TriggRate[i] = sum(TriggerAllSims[i])/len(TriggerAllSims[i])

for j in range(len(ZenithSortCut)):
    MeanTriggRate[j] = np.mean(TriggRate[ZenithAll== ZenithSortCut[j]])
    StdTriggRate[j] =  np.std(TriggRate[ZenithAll== ZenithSortCut[j]])
    
        
plt.scatter(np.array(ZenithAll)[TriggRate<1], TriggRate[TriggRate<1])    


plt.hist(ZenithAll, density=True)
plt.xlabel("zenith [Deg.]")
plt.savefig("/Users/chiche/Desktop/TZenithDistrib.pdf", bbox_inches = "tight")
plt.show()

plt.hist(EnergyAll, density=True, bins =100)
plt.xlabel("Energy [EeV]")
plt.savefig("/Users/chiche/Desktop/EDistrib.pdf", bbox_inches = "tight")
plt.show()


plt.errorbar(ZenithSortCut[1:], MeanTriggRate[1:], yerr= StdTriggRate[1:], fmt ='o', linestyle='-')
plt.xlabel('Zenith [Deg.]')
plt.ylabel("Trigger rate")
plt.title(r"$E = 10^{16.5} eV$") 
plt.grid()
plt.savefig("/Users/chiche/Desktop/TriggRate.pdf", bbox_inches = "tight")
plt.show()   
    
'''