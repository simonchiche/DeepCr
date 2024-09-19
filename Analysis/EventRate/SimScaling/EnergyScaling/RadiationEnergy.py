#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 22 02:21:35 2024

@author: chiche
"""


import numpy as np
import os
import subprocess
import matplotlib.pyplot as plt
import glob
import sys
import pickle
from scipy.optimize import curve_fit

sys.path.append("/Users/chiche/Desktop/DeepCrSearch/Analysis/")
from Modules.Fluence.FunctionsRadiationEnergy import GetRadiationEnergy, GetFluence
from Modules.Fluence.FunctionsGetFluence import GetDepths, CorrectLength, Traces_cgs_to_si
from datetime import datetime


font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)

Save = False

Energies = [0.0316, 0.1]
Zeniths = [20, 28, 34, 39, 43, 47, 50]
RadEnergyAll = dict()

#for i in range(Energies):
#    for j in range(Zeniths):
        
SimDir = "DeepCrLib"
SimName = "Rectangle_Proton_0.1_34_0_1"

# We create a directory where the outputs are stored
date = datetime.today().strftime('%Y-%m-%d')
WorkPath = os.getcwd()
OutputPath = WorkPath + "/Plots/" + SimDir + "/" + date + "/" 
cmd = "mkdir -p " + OutputPath
p =subprocess.Popen(cmd, cwd=os.getcwd(), shell=True)
stdout, stderr = p.communicate()
        

# =============================================================================
#                              Load Traces
# =============================================================================

Path =  "/Users/chiche/Desktop/DeepCrSearch"\
+ "/Simulations/" + SimDir + "/" + SimName + "/"
energy = float(Path.split("/")[-2].split("_")[2])
theta  = float(Path.split("/")[-2].split("_")[3])

SimDataPath = "/Users/chiche/Desktop/DeepCrSearch/Analysis/SimData/" + SimName 

if(not(os.path.exists(SimDataPath))):
    
    cmd = "mkdir -p " + SimDataPath
    p =subprocess.Popen(cmd, cwd=os.getcwd(), shell=True)
    stdout, stderr = p.communicate()
    Nant, Traces_C, Traces_G, Pos = LoadTraces(Path)
    
    np.save(SimDataPath + "/Nant", Nant)
    with open(SimDataPath + '/Traces_C.pkl', 'wb') as file:
        pickle.dump(Traces_C, file)
    with open(SimDataPath + '/Traces_G.pkl', 'wb') as file:
        pickle.dump(Traces_G, file)
    np.save(SimDataPath + "/Pos", Pos)
    
else:
    
    Nant  = np.load(SimDataPath + "/Nant.npy")
    with open(SimDataPath + '/Traces_C.pkl', 'rb') as file:
        Traces_C = pickle.load(file)
    with open(SimDataPath + '/Traces_G.pkl', 'rb') as file:
        Traces_G = pickle.load(file)
    Pos  = np.load(SimDataPath + "/Pos.npy", allow_pickle=True)
        
    
Nlay, Nplane, Depths = GetDepths(Pos)

# To resize the length of the Coreas traces
CorrectLength(Traces_C, False)

#cgs to si
Traces_C = Traces_cgs_to_si(Traces_C)
Traces_G = Traces_cgs_to_si(Traces_G)
    

# =============================================================================
#                         GetFluence
# =============================================================================

fx, fy, fz = GetFluence(Traces_G, Nant)
EradAll, Erad = GetRadiationEnergy(Traces_C, Nant, Nplane, Pos, Depths)


#15529314.152721744

#10^16.5

# 20: 1.59 /24
# 28: 2.8/ 32

#10^17
# 10: 8.52 / 78
# 20: 18 /104.24
# 28: 35 /148
# 34: 158

Erad16_5 = np.array([1.59, 2.8, 4.57])
Erad17 = np.array([18, 35, 40])
ScaledE17 = Erad16_5 *(1e17/10**16.5)**2

x16 = np.ones(len(Erad16_5))*16.5
x17 = np.ones(len(Erad17))*17

xscale = np.array([16.5, 17])

plt.scatter(x16, Erad16_5)
plt.scatter(x17, Erad17)
for i in range(len(x17)):
    yscale = np.array([Erad16_5[i], ScaledE17[i]])
    plt.plot(xscale, yscale)

plt.xlabel("Log10(E [eV])")
plt.ylabel("$E_{rad}$")
plt.title("In-air energy scaling", fontsize =12)
plt.show()




    
    