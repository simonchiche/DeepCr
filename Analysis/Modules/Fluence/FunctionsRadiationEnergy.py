#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep  5 15:51:27 2024

@author: chiche
"""


import numpy as np
import sys

sys.path.append("/Users/chiche/Desktop/DeepCrSearch/Analysis/")
import Modules.Fluence.FunctionsGetFluence


#from FunctionsGetFluence import Norm, LoadTraces, GetPeakTraces, \
#Traces_cgs_to_si, GetDepths, CorrectScaling, CombineTraces, CorrectLength, \
#GetIntTraces, GetIntTracesSum, GetRadioExtent


def GetFluence(TracesTot, Nant):
    
    
    eps0 = 8.85e-12
    c = 3e8
    
    fx, fy, fz = np.zeros(Nant), np.zeros(Nant), np.zeros(Nant)
    
    for i in range(Nant):
        if(i%10==0): print("%.d/%.d" %(i, Nant))
        
        DeltaT = TracesTot[i][-1, 0] - TracesTot[i][0, 0]
        #fx[i] = eps0*c*np.sum(TracesTot[i][:,1]**2)
        #fy[i] = eps0*c*np.sum(TracesTot[i][:, 2]**2)
        #fz[i] = eps0*c*np.sum(TracesTot[i][:, 3]**2)
        fx[i] = eps0*c*np.trapz(TracesTot[i][:,1]**2, TracesTot[i][:,0])
        fy[i] = eps0*c*np.trapz(TracesTot[i][:, 2]**2, TracesTot[i][:,0])
        fz[i] = eps0*c*np.trapz(TracesTot[i][:, 3]**2, TracesTot[i][:,0])
          
    return fx, fy, fz


def GetTimeFluence(TracesTot, Nant):
    
    
    eps0 = 8.85e-12
    c = 3e8
    
    f = np.zeros(Nant)
    
    for i in range(Nant):
        if(i%10==0): print("%.d/%.d" %(i, Nant))
        
        DeltaT = TracesTot[i][-1, 0] - TracesTot[i][0, 0]
        #fx[i] = eps0*c*np.sum(TracesTot[i][:,1]**2)
        #fy[i] = eps0*c*np.sum(TracesTot[i][:, 2]**2)
        #fz[i] = eps0*c*np.sum(TracesTot[i][:, 3]**2)
        I = eps0*c*TracesTot[i][:,1]**2 +\
                                TracesTot[i][:,2]**2 + eps0*c*TracesTot[i][:,3]**2
                                
        f[i] = np.trapz(I, TracesTot[i][:,0])
                                
    return np.array(f)



def GetRadiationEnergy(TracesTot, Nant, Nplane, Pos, Depths):
    
    fx, fy, fz = GetFluence(TracesTot, Nant)
   
    f =  GetTimeFluence(TracesTot, Nant)
    
    EradAll = np.zeros(len(Depths))
    
    for i in range(len(Depths)):
        
        step = abs(Pos[i*Nplane+ 1,0] - Pos[i*Nplane,0])
            
        for j in range(Nplane):
            
            index = j + i*Nplane
            #EradAll[i] = EradAll[i] + np.sqrt(fx[index] + fy[index] + fz[index])*step**2
            EradAll[i] =  EradAll[i] + f[index]*step**2
            
            if(index%100==0): print(step, Pos[index, 2])
    return EradAll, np.mean(EradAll)