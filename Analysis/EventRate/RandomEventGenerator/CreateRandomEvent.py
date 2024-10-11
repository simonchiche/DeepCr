#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Sep 20 12:12:01 2024

@author: chiche
"""
import os
import sys
import shutil
import random
import pickle
import numpy as np
import matplotlib.pyplot as plt
from SelectClosestEvent import LoadClosestEvent
sys.path.append("/Users/chiche/Desktop/DeepCrSearch/Analysis/")
from Modules.Shower.GetAntenaInfos import GetSurfaceAntennas, GetDepths, GetAntennaLayer

savedir = os.getcwd() + "/Plots/"

def GetRandomCore(AntPos, glevel, Nrand):
    
    SurfaceAntennas =  GetSurfaceAntennas(AntPos, glevel)
    
    xmin, xmax =  min(SurfaceAntennas[:,0]), max(SurfaceAntennas[:,0])
    ymin, ymax =  min(SurfaceAntennas[:,1]), max(SurfaceAntennas[:,1])
    
    xrand = np.random.uniform(xmin, xmax, Nrand)
    yrand = np.random.uniform(ymin, ymax, Nrand)
    
    return xrand, yrand, SurfaceAntennas

def GetRandomCoreCircular(rlim, Nrand):
    
    rmin = 0
    phimin = 0
    phimax = 360
    rrand = np.sqrt(np.random.uniform(rmin, rlim**2, Nrand))
    phirand = np.random.uniform(phimin, phimax, Nrand)
    xrand = rrand*np.cos(phirand*np.pi/180.0)
    yrand = rrand*np.sin(phirand*np.pi/180.0)
   
    return xrand, yrand, rrand

def GetDeepLayerCore(AntPos,glevel, Depth):
    
    DeepAntennas =  GetAntennaLayer(AntPos, glevel, Depth)
    xcore,ycore = DeepAntennas[:, 0], DeepAntennas[:, 1]
    
    return xcore, ycore
    

def ShiftAntennaPos(AntPos, xrand, yrand, phi_rand, glevel):
    
    # Rotation of the antenna positions
    rotation_matrix = np.array([[np.cos(phi_rand), -np.sin(phi_rand)], 
                                [np.sin(phi_rand), np.cos(phi_rand)]])
    
    xy_pos = np.copy(AntPos[:, :2])
    rotated_xy = np.dot(xy_pos, rotation_matrix.T)
    
    # Shifting the antennas to the core impact postion
    rotated_xy[:, 0] += xrand 
    rotated_xy[:, 1] += yrand
    
    plt.scatter(rotated_xy[:, 0], rotated_xy[:, 1])
    plt.scatter(xrand, yrand)
    plt.show()
    #sys.exit()
    ShiftedAntPos = np.copy(AntPos)
    ShiftedAntPos[:, :2] = rotated_xy     

    return ShiftedAntPos

def RotateAntennaPos(AntPos, phi_rand):
    
    # Rotation of the antenna positions
    rotation_matrix = np.array([[np.cos(phi_rand), -np.sin(phi_rand)], 
                                [np.sin(phi_rand), np.cos(phi_rand)]])
    
    xy_pos = np.copy(AntPos[:, :2])
    rotated_xy = np.dot(xy_pos, rotation_matrix.T)

    #plt.scatter(rotated_xy[:, 0], rotated_xy[:, 1])
    #plt.show()
    #sys.exit()
    RotatedAntPos = np.copy(AntPos)
    RotatedAntPos[:, :2] = rotated_xy     

    return RotatedAntPos

def GetClosestPos(AntPos, xtarget, ytarget):
    
    r_all = np.sqrt((AntPos[:,0] - xtarget)**2 + (AntPos[:,1] - ytarget)**2)
    i_sel = np.argmin(r_all)
    
    x_sel, y_sel = AntPos[i_sel,0],  AntPos[i_sel,1]
    
    return x_sel, y_sel, i_sel
    
    
def SelectPowerLineEvent(E_rand, theta_rand, phi_rand, glevel, SimPath, DataPath):
    
    Traces_C, Traces_G, AntPos, Selfile = \
    LoadClosestEvent(E_rand, theta_rand, phi_rand, SimPath, DataPath)
    
    Nlay, NantLay, Depths = GetDepths(AntPos, glevel)
    
    xcore, ycore,  = \
    GetRandomCore(AntPos, glevel, 1)[0:2]

    # Shifted antenna positions
    ShiftedPos = ShiftAntennaPos(AntPos, xcore, ycore, phi_rand, glevel)
    Plot = True
    if(Plot):
        plt.scatter(ShiftedPos[:729,0], ShiftedPos[:729,1], s= 1)
        plt.scatter(xcore, ycore)
    
    print(Nlay, NantLay, Depths)
    
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
        
    return AntPow, Traces_C_pow, Traces_G_pow, xcore, ycore



def GetDeepAntennaLayerEvents(E_rand, theta_rand, phi_rand, glevel, SimPath, DataPath):
    
    Traces_C, Traces_G, AntPos, Selfile = \
    LoadClosestEvent(E_rand, theta_rand, phi_rand, SimPath, DataPath)
    
    Nlay, NantLay, Depths = GetDepths(AntPos, glevel)
    #print(Nlay, NantLay, Depths)
    
    # Shifted antenna positions
    RotatedPos = RotateAntennaPos(AntPos, phi_rand)
    
    xdeep, ydeep,  = \
    GetDeepLayerCore(RotatedPos, glevel, Depths[0])
    
    Plot = False
    if(Plot):
        plt.scatter(RotatedPos[:729,0], RotatedPos[:729,1], s= 1)
        plt.scatter(xdeep, ydeep)
    
    indices = np.where(RotatedPos[:,2] == (glevel- Depths[0]))[0]
    
    Traces_C_pow = {i: Traces_C[i] for i in indices}
    Traces_G_pow = {i: Traces_G[i] for i in indices}
    
    zdeep = AntPos[indices,2]
        
    AntPow = np.array([xdeep, ydeep, zdeep]).T
        
    return AntPow, Traces_C_pow, Traces_G_pow, xdeep, ydeep, Selfile

def generate_antenna_file(x_positions, y_positions, z_positions, path, filename="SIM000001.list"):
    # check array size
    if not (len(x_positions) == len(y_positions) == len(z_positions)):
        raise ValueError("Les listes x, y, et z doivent avoir la même longueur.")
    
    full_filename = path +"/" + filename
    with open(full_filename, 'w') as file:
        for i, (x, y, z) in enumerate(zip(x_positions, y_positions, z_positions)):
            file.write(f"AntennaPosition = {x:.1f} {y:.1f} {z:.1f} ch{i}\n")

    return

def save_traces(Traces):
    """
    Save each numpy array from the Traces dictionary into a separate file.
    
    Args:
    Traces (dict): A dictionary where each key corresponds to a numpy array [N, 3].
    """
    
    # Iterate through the dictionary keys and save each array to a file
    for i in Traces:
        filename = f"raw_ch{i}.dat"
        np.savetxt(filename, Traces[i], fmt="%.8e\t%.8e\t%.8e")
        print(f"Saved {filename}")
        
    return
    
def SaveEvent(AntPos, Traces_C, Traces_G, dirname, LibPath, j):
    
    savedir = LibPath + "/" + dirname + "_" + str(j)
    print(savedir)
    try:
        os.mkdir(savedir)
    except FileExistsError:
        print("Directory already exists.")
    os.chdir(savedir)
    CoreasFolder = savedir + "/" + "SIM000001_coreas"
    GeantFolder = savedir + "/" + "SIM000001_geant"
    try:
        os.mkdir(CoreasFolder)
        os.mkdir(GeantFolder)
    except FileExistsError:
        print("Directory already exists.")
    
    os.chdir(CoreasFolder)
    save_traces(Traces_C)
    
    os.chdir(GeantFolder)
    save_traces(Traces_G)
        
    #with open(CoreasFolder + '/Traces_C.pkl', 'wb') as file:
    #    pickle.dump(Traces_C, file)
    #with open(GeantFolder + '/Traces_G.pkl', 'wb') as file:
    #    pickle.dump(Traces_G, file)
            
    #np.save(savedir + "/Pos", AntPos)
    generate_antenna_file(AntPos[:,0], AntPos[:,1], AntPos[:,2], \
                          savedir, filename="SIM000001.list")
    
    source = "/Users/chiche/Desktop/SimFiles/" + str(dirname.split("_")[2]) + \
    "/" + str(dirname.split("_")[3])  + "/"
    print(source)

    for file_name in os.listdir(source):
        full_file_name = os.path.join(source, file_name)
        if os.path.isfile(full_file_name):
            shutil.copy(full_file_name, savedir)
    return       
    
Test = False
if(Test):

    SimPath = "/Users/chiche/Desktop/DeepCrSearch/Simulations/DeepCrLib"
    DataPath = "/Users/chiche/Desktop/DeepCrSearch/Analysis/SimData"
    glevel =3216
    AntPos = \
    LoadClosestEvent(0.1, 10, 24, SimPath, DataPath)[2]
    xrand, yrand, SurfaceAntennas =\
    GetRandomCore(AntPos, glevel, 1)
    
    xrandcirc, yrandcirc, rrand = GetRandomCoreCircular(500, 10)
    
    plt.scatter(SurfaceAntennas[:,0], SurfaceAntennas[:,1], label ="antennas", s= 15)
    plt.scatter(xrandcirc, yrandcirc, label ="random core", s= 15)
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.legend()
    if(True): plt.savefig(savedir +"RandomCoresCirc2.pdf", bbox_inches = "tight")
    plt.show()
    
    xrand, yrand,  = \
    GetRandomCore(AntPos, glevel, 1)[0:2]
    
    ShiftedPos = \
    ShiftAntennaPos(AntPos, xrand, yrand, 24, glevel)
    
    plt.scatter(ShiftedPos[:729,0], ShiftedPos[:729,1], s=1, \
                label =r"shifted pos ($\varphi_{\rm rand} = 45^{\circ})$")
    plt.scatter(xrand, yrand, label ="random core")
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    plt.legend()
    if(False): plt.savefig(savedir +"ShiftedArray.pdf", bbox_inches = "tight")
    plt.show()
    
    AntPow, Traces_C_pow, Traces_G_pow, xrand, yrand = \
    SelectPowerLineEvent(0.1, 10, 24, glevel, SimPath, DataPath)
    
    AntPow, Traces_C_pow, Traces_G_pow, xrand, yrand = \
    GetDeepAntennaLayerEvents(0.1, 10, 24, glevel, SimPath, DataPath)
    
    plt.plot(AntPow[:,0], AntPow[:,2], label ="power line vs x-axis")
    plt.scatter(AntPow[:,0], AntPow[:,2])
    plt.scatter(xrand, glevel, color = 'red', label ="core")
    plt.xlabel("x [m]")
    plt.ylabel("z [m]")
    plt.legend()
    if(False): plt.savefig(savedir +"PoweLinex.pdf", bbox_inches = "tight")
    plt.show()
    
    plt.plot(AntPow[:,1], AntPow[:,2], label ="power line vs y-axis")
    plt.scatter(AntPow[:,1], AntPow[:,2])
    plt.scatter(yrand, glevel, color = 'red', label ="core")
    plt.xlabel("y [m]")
    plt.ylabel("z [m]")
    plt.legend()
    if(False): plt.savefig(savedir +"PoweLiney.pdf", bbox_inches = "tight")
    plt.show()
    
    plt.scatter(AntPow[:,0], AntPow[:,1])
    plt.xlabel("x [m]")
    plt.ylabel("y [m]")
    if(False): plt.savefig(savedir +"PoweLinexy.pdf", bbox_inches = "tight")
    plt.show()
    