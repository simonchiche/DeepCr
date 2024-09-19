#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 23:15:45 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt
import os
import sys

savedir = os.getcwd() + "/Plots/"


# We want to generate a flat distribution in cos(theta) 
#between 0 and 50 degrees and a uniform distribution in azimuth angle


def GetRandomTheta(min_th, max_th, N):
    
    min_cos_theta = np.cos(min_th*180/np.pi)
    max_cos_theta = np.cos(max_th*np.pi/180.0)
    cos_theta_all = np.random.uniform(min_cos_theta, max_cos_theta, N)
    theta_all = np.arccos(cos_theta_all)*180/np.pi

    return np.round(theta_all, 1)


def GetRandomPhi(min_phi, max_phi, N):
    
    phi_all = np.random.uniform(min_phi, max_phi, N)

    return np.round(phi_all,1)


# =============================================================================
#                               Test plots
# =============================================================================

Plot = False
if(Plot):
    
    min_th = 0
    max_th = 50 
    N = 5000
    theta_all = GetRandomTheta(min_th, max_th, N)
    plt.hist(theta_all)
    plt.xlabel("zenith [Deg.]")
    plt.ylabel("counts")
    plt.savefig(savedir + "theta_distrib.pdf", bbox_inches = "tight")
    plt.show()
    
    
    min_phi = 0
    max_phi = 360
    phi_all = GetRandomPhi(min_phi, max_phi, N)
    plt.hist(phi_all)
    plt.xlabel("azimuth [Deg.]")
    plt.ylabel("counts")
    plt.savefig(savedir + "phi_distrib.pdf", bbox_inches = "tight")
    plt.show()