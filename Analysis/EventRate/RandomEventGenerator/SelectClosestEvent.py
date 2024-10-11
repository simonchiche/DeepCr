#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 11:00:50 2024

@author: chiche
"""

import sys
import numpy as np
import glob
import pickle
import matplotlib.pyplot as plt

def GetSimParameters(SimPath):
    
    Simfiles =  glob.glob(SimPath + "/*")
    Nsim = len(Simfiles)
    E_sim = np.zeros(Nsim)
    theta_sim = np.zeros(Nsim)
    phi_sim = np.zeros(Nsim)
    
    
    for i in range(Nsim):
        
        splitfile = Simfiles[i].split('_')
        E_sim[i] = splitfile[2]
        theta_sim[i] = splitfile[3]
        phi_sim[i] = splitfile[4] 
    
    return E_sim, theta_sim, phi_sim, Simfiles


def GetClosestEvent(E_rand, theta_rand, phi_rand, SimPath):
    
    E_sim, theta_sim, phi_sim, Simfiles = GetSimParameters(SimPath)
    
    E_sel = E_sim[np.argmin(abs(E_sim-E_rand))]
    theta_sel = theta_sim[np.argmin(abs(theta_sim - theta_rand))]
    phi_sel = phi_sim[np.argmin(abs(phi_sim - phi_rand))]
    SelParams = np.array([E_sel, theta_sel, phi_sel])
    
    splitfile = Simfiles[0].split("/")[-1].split('_')
    splitfile[2] = str(SelParams[0])
    splitfile[3] = str(int(SelParams[1]))
    splitfile[4] = str(int(SelParams[2]))
    
    Selfile = '_'.join(splitfile)
        
    return Selfile


def LoadClosestEvent(E_rand, theta_rand, phi_rand, SimPath, DataPath):
    
    Selfile = GetClosestEvent(E_rand, theta_rand, phi_rand, SimPath)
    print(Selfile)
    SimDataPath = DataPath + "/" + Selfile 
    
    with open(SimDataPath + '/Traces_C.pkl', 'rb') as file:
        Traces_C = pickle.load(file)
    with open(SimDataPath + '/Traces_G.pkl', 'rb') as file:
        Traces_G = pickle.load(file)
    AntPos  = np.load(SimDataPath + "/Pos.npy", allow_pickle=True)
    
    return Traces_C, Traces_G, AntPos, Selfile


# Testing of the code
Test = False
if(Test):
    SimPath = "/Users/chiche/Desktop/DeepCrSearch/Simulations/DeepCrLib"
    data = GetSimParameters(SimPath)
    Selfile = GetClosestEvent(0.01, 46, 34, SimPath)
    DataPath = "/Users/chiche/Desktop/DeepCrSearch/Analysis/SimData"
    
    Traces_C, Traces_G, AntPos = \
    LoadClosestEvent(0.01, 46, 34, SimPath, DataPath)

