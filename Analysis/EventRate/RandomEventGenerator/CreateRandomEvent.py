#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 12:12:01 2024

@author: chiche
"""
import os
import sys
import random
import numpy as np
import matplotlib.pyplot as plt
from SelectClosestEvent import LoadClosestEvent
from Modules.Shower.GetAntenaInfos import GetSurfaceAntennas, GetDepths, GetAntennaLayer

savedir = os.getcwd() + "/Plots/"

def GetRandomCore(AntPos, glevel, Nrand):
    
    SurfaceAntennas =  GetSurfaceAntennas(AntPos, glevel)
    
    xmin, xmax =  min(SurfaceAntennas[:,0]), max(SurfaceAntennas[:,0])
    ymin, ymax =  min(SurfaceAntennas[:,1]), max(SurfaceAntennas[:,1])
    
    xrand = np.random.uniform(xmin, xmax, Nrand)
    yrand = np.random.uniform(ymin, ymax, Nrand)
    
    return xrand, yrand, SurfaceAntennas


def ShiftAntennaPos(AntPos, phi_rand, glevel):
    
    xrand, yrand,  = \
    GetRandomCore(AntPos, glevel, 1)[0:2]
    
    # Rotation of the antenna positions
    rotation_matrix = np.array([[np.cos(phi_rand), -np.sin(phi_rand)], 
                                [np.sin(phi_rand), np.cos(phi_rand)]])
    
    xy_pos = np.copy(AntPos[:, :2])
    rotated_xy = np.dot(xy_pos, rotation_matrix.T)
    
    # Shifting the antennas to the core impact postion
    rotated_xy[:, 0] += xrand 
    rotated_xy[:, 1] += yrand
    
    ShiftedAntPos = np.copy(AntPos)
    ShiftedAntPos[:, :2] = rotated_xy     

    return ShiftedAntPos, xrand, yrand

def GetClosestPos(AntPos, xtarget, ytarget):
    
    r_all = np.sqrt((AntPos[:,0] - xtarget)**2 + (AntPos[:,1] - ytarget)**2)
    i_sel = np.argmin(r_all)
    
    x_sel, y_sel = AntPos[i_sel,0],  AntPos[i_sel,1]
    
    return x_sel, y_sel, i_sel
    
    
def SelectPowerLineEvent(E_rand, theta_rand, phi_rand, glevel, SimPath, DataPath):
    
    Traces_C, Traces_G, AntPos = \
    LoadClosestEvent(E_rand, theta_rand, phi_rand, SimPath, DataPath)
    
    # Shifted antenna positions
    ShiftedPos, xrand, yrand = ShiftAntennaPos(AntPos, phi_rand, glevel)
    
    Plot = True
    if(Plot):
        plt.scatter(ShiftedPos[:729,0], ShiftedPos[:729,1], s= 1)
        plt.scatter(xrand, yrand)
    
    
    Nlay, NantLay, Depths = GetDepths(ShiftedPos, glevel)

    
    # We aim to find the antennas that are the closest 
    #to the power line position in (0,0)
    xtarget, ytarget = 0, 0 
    
    Traces_C_pow = dict()
    Traces_G_pow = dict()
    xpow, ypow, zpow = np.zeros(Nlay), np.zeros(Nlay), np.zeros(Nlay)
    
    for i in range(Nlay):
        
        # For each layer we get the (x,y) coordinated of the antenna positions
        xy_lay = GetAntennaLayer(ShiftedPos, glevel, Depths[i])[:,:2]
        # We select the antennas that are the closest to the power line
        xsel, ysel, isel = GetClosestPos(xy_lay, xtarget, ytarget)
        xpow[i], ypow[i] = xsel, ysel
        zpow[i] = AntPos[i*NantLay,2]
        isel_lay = isel + i*NantLay
        Traces_C_pow[i] = Traces_C[isel_lay]
        Traces_G_pow[i] = Traces_G[isel_lay]
        
        # Since our antennas are not aligned along the z-direction, 
        # we add this condition to get  a string of antenna which are the
        # closest to the one selected at a 100 m depth
        if(i==0): xtarget, ytarget = xsel, ysel
    if(Plot):
        plt.scatter(xtarget, ytarget)
        plt.show()
        
    AntPow = np.array([xpow, ypow, zpow]).T
        
    return AntPow, Traces_C_pow, Traces_G_pow, xrand, yrand
          
    
Test = False
if(Test):

    SimPath = "/Users/chiche/Desktop/DeepCrSearch/Simulations/DeepCrLib"
    DataPath = "/Users/chiche/Desktop/DeepCrSearch/Analysis/SimData"
    glevel =3216
    AntPos = \
    LoadClosestEvent(0.1, 10, 24, SimPath, DataPath)[2]
    xrand, yrand, SurfaceAntennas =\
    GetRandomCore(AntPos, glevel, 1)
    
    plt.scatter(SurfaceAntennas[:,0], SurfaceAntennas[:,1], label ="antennas", s= 15)
    plt.scatter(xrand, yrand, label ="random core", s= 15)
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.legend()
    if(True): plt.savefig(savedir +"RandomCores.pdf", bbox_inches = "tight")
    plt.show()
    
    ShiftedPos, xrand, yrand = \
    ShiftAntennaPos(AntPos, 24, glevel)
    
    plt.scatter(ShiftedPos[:729,0], ShiftedPos[:729,1], s=1, \
                label =r"shifted pos ($\varphi_{\rm rand} = 45^{\circ})$")
    plt.scatter(xrand, yrand, label ="random core")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.legend()
    if(True): plt.savefig(savedir +"ShiftedArray.pdf", bbox_inches = "tight")
    plt.show()
    
    AntPow, Traces_C_pow, Traces_G_pow, xrand, yrand = \
    SelectPowerLineEvent(0.1, 10, 24, glevel, SimPath, DataPath)
    
    plt.plot(AntPow[:,0], AntPow[:,2], label ="power line vs x-axis")
    plt.scatter(AntPow[:,0], AntPow[:,2])
    plt.scatter(xrand, glevel, color = 'red', label ="core")
    plt.xlabel("x [m]")
    plt.ylabel("z [m]")
    plt.legend()
    if(True): plt.savefig(savedir +"PoweLinex.pdf", bbox_inches = "tight")
    plt.show()
    
    plt.plot(AntPow[:,1], AntPow[:,2], label ="power line vs y-axis")
    plt.scatter(AntPow[:,1], AntPow[:,2])
    plt.scatter(yrand, glevel, color = 'red', label ="core")
    plt.xlabel("y [m]")
    plt.ylabel("z [m]")
    plt.legend()
    if(True): plt.savefig(savedir +"PoweLiney.pdf", bbox_inches = "tight")
    plt.show()
    
    plt.scatter(AntPow[:,0], AntPow[:,1])
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    if(True): plt.savefig(savedir +"PoweLinexy.pdf", bbox_inches = "tight")
    plt.show()