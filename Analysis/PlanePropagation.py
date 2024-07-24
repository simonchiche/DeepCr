#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 03:46:27 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt
from IceCorePosition import core
from GetXmaxCoordinate import getXmaxPosition
from IncidentAngles import GetIncidentAngles
from FunctionsPlotFluence import EfieldMap
from FunctionsGetFluence import GetIntTraces, GetDepths
import sys
import os
import pickle
from datetime import datetime


Xmax = 765.70844
phi = 0
theta = 70
glevel =  3216
injh =1e5
depth = 100
model = 1
SimDir = "InclinedShower"
SimName = "Rectangle_Proton_1.0_70_0_1"
energy = float(SimName.split("_")[2])
SimDataPath = "/Users/chiche/Desktop/DeepCrSearch/Analysis/SimData/" + SimName 
AntPath = SimDataPath + "/Pos.npy"

date = datetime.today().strftime('%Y-%m-%d')
WorkPath = os.getcwd()
OutputPath = WorkPath + "/Plots/" + SimDir + "/" + date + "/" 

#  Loading antennas
AntPos  = np.load(AntPath, allow_pickle=True)
Nlay, Nplane, Depths = GetDepths(AntPos)
SurfacePos = np.copy(AntPos)[Nplane:]

#  Loading Efield
with open(SimDataPath + '/Traces_tot.pkl', 'rb') as file:
    Traces_tot = pickle.load(file)    



def GetCurvedDist(xarr, zarr):
    
    xarr, zarr = np.array(xarr), np.array(zarr)
    
    xarr = xarr[zarr<=0]
    zarr = zarr[zarr<=0]
    
    points = np.array([xarr, zarr]).T
    differences = np.diff(points, axis=0)
    distances = np.linalg.norm(differences, axis=1)

    CurvedDist = np.sum(distances)
    
    return CurvedDist

def PlanePropagation(Xmax, phi, theta, glevel, SurfacePos, depth, Traces_tot, model):
    
    XmaxPos = getXmaxPosition(Xmax, phi, theta, glevel)
    
    theta_i, phi_i  = GetIncidentAngles(XmaxPos, SurfacePos)

    for i in range(len(theta_i)):
        
        x_ray, z_ray = core(depth, theta_i[i], model)[0:2]
        rcore = x_ray[-1]
        dpath = GetCurvedDist(x_ray, z_ray)
        #print(dpath, rcore)
        
        SurfacePos[i][0] = SurfacePos[i][0] + rcore*np.cos(phi*np.pi/180.0)
        SurfacePos[i][1] = SurfacePos[i][1] + rcore*np.sin(phi*np.pi/180.0)
        Traces_tot[i][:,1:] = Traces_tot[i][:,1:]/dpath
        
    
    SurfacePos[:, 2] = SurfacePos[:, 2] - depth
    DepthPos = SurfacePos
    
    return np.array(DepthPos), Traces_tot
  
DepthPos, Traces_tot = PlanePropagation(Xmax, phi, theta, glevel, SurfacePos, depth,Traces_tot, model)  

Ex_tot_int, Ey_tot_int, Ez_tot_int, Etot_int = GetIntTraces(Traces_tot, len(Traces_tot))

Plot = True
if(Plot):
    plt.scatter(AntPos[Nplane:][:, 0], AntPos[Nplane:][:, 1])    
    plt.scatter(DepthPos[:, 0], DepthPos[:, 1])  
    plt.show()
    
    EfieldMap(DepthPos, 1, Nplane, Etot_int[Nplane:], \
              "InclinedScaled", True, energy, theta, OutputPath)
    
    EfieldMap(DepthPos, 1, Nplane, Etot_int[Nplane:], \
          "InclinedScaled", False, energy, theta, OutputPath)

