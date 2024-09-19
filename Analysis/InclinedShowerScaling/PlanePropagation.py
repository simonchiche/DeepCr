#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 24 03:46:27 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt
import sys
sys.path.append("/Users/chiche/Desktop/DeepCrSearch/Analysis/")
from Modules.RayTracing.IceCorePosition import core
from Modules.Shower.GetXmaxCoordinate import getXmaxPosition 
from Modules.Shower.IncidentAngles import GetIncidentAngles 
from Modules.Fluence.FunctionsPlotFluence import EfieldMap
from Modules.Fluence.FunctionsGetFluence import GetIntTraces, GetDepths
import os
import pickle
from datetime import datetime
import subprocess

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
cmd = "mkdir -p " + OutputPath
p =subprocess.Popen(cmd, cwd=os.getcwd(), shell=True)
stdout, stderr = p.communicate()

#  Loading antennas
AntPos  = np.load(AntPath, allow_pickle=True)
Nlay, Nplane, Depths = GetDepths(AntPos)
SurfacePos = np.copy(AntPos)[Nplane:]

#  Loading Efield
with open(SimDataPath + '/Traces_tot.pkl', 'rb') as file:
    Traces_tot = pickle.load(file)    

Ex_tot_sim, Ey_tot_sim, Ez_tot_sim, Etot_sim = GetIntTraces(Traces_tot, len(Traces_tot))


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



# =============================================================================
#                       Comparison with simulations
# =============================================================================



condition = (AntPos[:Nplane, 0] >  -2713.3) & (AntPos[:Nplane, 0] <  2995.3) \
& (AntPos[:Nplane, 1] >  -2950.7) & (AntPos[:Nplane, 1] <  2850.7)

PosCut =np.array([AntPos[:Nplane, 0][condition], AntPos[:Nplane, 1][condition], \
                  AntPos[:Nplane, 2][condition]]).T

Etot_simcut = Etot_sim[:Nplane][condition]
    
plt.scatter(PosCut[:,0], PosCut[:,1], c= Etot_simcut, cmap ="jet")
cbar = plt.colorbar()
plt.xlabel("x [m]")
plt.ylabel("y [m]")
cbar.set_label("E [$\mu V/m$]")
plt.legend(["Depth = %.f m" %(100)], loc ="upper right")
plt.title("SimCut" + " map (E =$%.2f\,$EeV, $\\theta=%.1f^{\circ}$)" %(energy, theta), size =14)
#plt.savefig\
#            (OutputPath + "SimCut" + "EfieldMap_E%.2f_th%.1f_depth%1.f.pdf" \
#             %(energy, theta, 100), bbox_inches = "tight")
plt.show()

# Comparison of the traces at the antenna level

Error = np.zeros(len(PosCut))

for i in range(len(PosCut)):
    
    points_sim = np.array([PosCut[i,0], PosCut[i,1]]).T
    
    points_scaled = np.array([DepthPos[:,0], DepthPos[:,1]]).T
    
    distances = np.linalg.norm(points_scaled - points_sim, axis=1) 
    
    kmin = np.argmin(distances)
    
    Error[i] = (Etot_simcut[i] - Etot_int[Nplane:][kmin])/Etot_simcut[i]
    

plt.plot(Error)
plt.xlabel("Antenna ID")
plt.ylabel("(Sim-Scaled)/Sim")
#plt.savefig\
#            (OutputPath + "ScaledInclined" + "RelErr_E%.2f_th%.1f_depth%1.f.pdf" \
#             %(energy, theta, 100), bbox_inches = "tight")
plt.show()    


plt.scatter(PosCut[:,0], PosCut[:,1], c= Error, cmap ="jet")
cbar = plt.colorbar()
plt.xlabel("x [m]")
plt.ylabel("y [m]")
cbar.set_label("(Sim-Scaled)/Sim")
plt.legend(["Depth = %.f m" %(100)], loc ="upper right")
plt.title("SimCut" + " map (E =$%.2f\,$EeV, $\\theta=%.1f^{\circ}$)" %(energy, theta), size =14)
#plt.savefig\
#            (OutputPath + "ScaledInclined" + "RelErrMap_E%.2f_th%.1f_depth%1.f.pdf" \
#             %(energy, theta, 100), bbox_inches = "tight")
plt.show()
    
    
plt.scatter(PosCut[:,0], PosCut[:,1], c= abs(Error), cmap ="jet")
cbar = plt.colorbar()
plt.xlabel("x [m]")
plt.ylabel("y [m]")
cbar.set_label("(Sim-Scaled)/Sim")
plt.legend(["Depth = %.f m" %(100)], loc ="upper right")
plt.title("SimCut" + " map (E =$%.2f\,$EeV, $\\theta=%.1f^{\circ}$)" %(energy, theta), size =14)
#plt.savefig\
#            (OutputPath + "ScaledInclined" + "RelErrAbsMap_E%.2f_th%.1f_depth%1.f.pdf" \
#             %(energy, theta, 100), bbox_inches = "tight")
plt.show()    
    

IntSim = 0

for i in range(len(Etot_simcut)):
    
    IntSim = IntSim + Etot_simcut[i]*1.0/len(Etot_simcut)
 
    
IntScaled = 0

for i in range(len(Etot_int)):
    
    IntScaled = IntScaled + Etot_int[i]*1.0/len(Etot_int)
    
    
    
    
    