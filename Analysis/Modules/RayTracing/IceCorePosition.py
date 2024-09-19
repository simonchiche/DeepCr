#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr  2 05:49:29 2024

@author: chiche
"""

####
# Code to estimate where to place the center of an antenna l

####
import numpy as np
import matplotlib.pyplot as plt
import sys

font = {'family' : 'normal',
        'weight' : 'normal',
        'size'   : 14}

plt.rc('font', **font)
theta = 30
depth = 100
dz = 1

def deg2rad(theta):
    return theta*np.pi/180.0

def n(z, model):
   
    depth = abs(z)
    if(model ==1): # Greenland
        A = 1.775
        if(z>-14.9):
            #print(z)
            B = -0.5019
            C =0.03247 
           
            n = A + B*np.exp(-C*depth)
        else:
            #print("here")
            B = -0.448023
            C = 0.02469 
            
            n = A + B*np.exp(-C*depth)
        
        if(z>0): n =1
        
    
    if(model == 2):
        #print("here")
        A = 1.775
        B = -0.43
        C =0.0132
        
        n = A + B*np.exp(-C*depth)
    if(z>0): n =1
        
    return n


def core(depth, theta, model):
   
    dz =1
    i1 = deg2rad(theta)
    theta = deg2rad(theta)
    height = 1e2
    x0  =  -height*np.tan(theta)
    z0 = height
    
    xall, zall, nall = [x0], [z0], [1]
    xlin, zlin = xall.copy(), zall.copy()
    
    z = height
        
    while(z>-depth):
    
        xold = xall[-1]
        zold = zall[-1]
        n1 = n(z, model)
        n2 = n(z - dz, model)
        #print(n1)
        #print(n1,n2, z)
        
        i2 = np.arcsin(n1/n2*np.sin(i1))
        i1 = i2

        #print(i2*180/np.pi)
        
        #print(i2*180/np.pi)
        xnew = xold + dz*np.tan(i2)
        znew = zold - dz
        
        xall.append(xnew)
        zall.append(znew)
        xlin.append(xlin[-1] + dz*np.tan(theta))
        zlin.append(zlin[-1] - dz)
        nall.append(n2)
        
        z = z -dz
    
    return xall, zall, xlin, zlin, nall


# =============================================================================
#                         Pllot test
# =============================================================================
    
Plot = False

if(Plot):

    xarr, zarr, xlin, zlin, narr = core(depth, theta, 2)
    
    # Greenland: curved 30º vs linear
    plt.plot(xarr, zarr, label = r"Greenland: $\theta = 30^{\circ}$")
    plt.plot(xlin, zlin, linestyle ="--", c="red", label ="linear")
    plt.xlabel("x [m]")
    plt.ylabel("z [m]")
    plt.legend()
    #plt.savefig("/Users/chiche/Desktop/GreenlandVsLin_theta30.pdf", bbox_inches='tight')
    plt.show()
    
    
    # Greenland: ray-bending  vs theta
    thetaAll = np.arange(0,60,10)
    for i in range(len(thetaAll)):
        
        xray, zray= core(depth, thetaAll[i], 1)[0:2]
        plt.plot(xray, zray, label = "$%.d^{\circ}$" %thetaAll[i])
        #plt.xlim(-50,200)
        plt.xlabel("x [m]")
        plt.ylabel("z [m]")
        plt.legend()
        plt.title("ray-bending Greenland", fontsize =14)
        #plt.xlim(-10,10)
    
    #plt.savefig("/Users/chiche/Desktop/RayBendingvsTheta.pdf", bbox_inches='tight')    
    plt.show()
    
    
    
    
    
    '''
    # SouthPole: ray-bending  vs theta
    for i in range(len(thetaAll)):
        
        xray, zray= core(depth, thetaAll[i], 2)[0:2]
        plt.plot(xray, zray, label = "$%.d^{\circ}$" %thetaAll[i])
        #plt.xlim(-50,200)
        plt.xlabel("x [m]")
        plt.ylabel("z [m]")
        plt.legend()
        plt.title("ray-bending South Pole", fontsize =14)
    
    #plt.savefig("/Users/chiche/Desktop/RayBendingvsThetaSpole.pdf", bbox_inches='tight')    
    plt.show()
    
    
    # Greenland - SouthPole: ray-bending  vs theta
    for i in range(len(thetaAll)):
        
        xray, zray= core(depth, thetaAll[i], 1)[0:2]
        xray_spole, zray_spole = core(depth, thetaAll[i], 2)[0:2]
        plt.plot(np.array(xray)-np.array(xray_spole),  np.array(zray), label = "$%.d^{\circ}$" %thetaAll[i])
        #plt.xlim(-50,200)
        plt.xlabel(r"$\Delta x $[m]")
        plt.ylabel("z [m]")
        plt.legend()
        plt.title("ray-bending Greenland - South Pole", fontsize =14)
    #plt.savefig("/Users/chiche/Desktop/RayBendingSpolevsGreenland.pdf", bbox_inches='tight')    
    plt.show()
    '''


###k(depth, theta, energy)
# check xmax, coreas files

