#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 17:28:26 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

### Script that generates random energies following a given energy spectrum

# Auger energy spectrum from PoS(ICRC2021)324
def EspectrumAuger(E):
    
    gamma_i = np.array([3.09, 2.85, 3.283, 2.54, 3.03, 5.3])   
    E_ij = np.array([2.8e16, 1.58e17, 5.0e18, 1.4e19, 4.7e19])  
    w_ij = [0.25,0.25, 0.05, 0.05, 0.05]  
    J0 = 8.34e-11
    
    J = J0*(E/1e16)**(-gamma_i[0])
    
    for i in range(len(E_ij)):
        
        J = J*(1+(E/E_ij[i])**(1.0/w_ij[i]))**\
        ((gamma_i[i] -  gamma_i[i+1])*w_ij[i])
        
    return J

# Generating a random energy following a spectrum
def generate_random_energy(CDF, E):
    # Random number between 0 and 1
    u = np.random.rand()
    
    # Finding the index corresping to u in the CDF
    energy_index = np.searchsorted(CDF, u) - 1
    
    # We return the corresponding energy
    return E[energy_index]

def GetRandomE(pmin, pmax, N):
    
    E = np.logspace(pmin, pmax, 1000)
    
    # Total integral for normalization
    Ntot = quad(EspectrumAuger, min(E), max(E))[0]

    # Cummulative Distribution Function
    CDF = np.zeros(len(E))
    
    for i in range(1, len(E)):
        # Cummulative integral from Emin to E[i]
        CDF[i] = quad(EspectrumAuger, min(E), E[i])[0] / Ntot

    # Generating random energies following the spectrum
    random_energies = [generate_random_energy(CDF, E) for _ in range(N)]
    
    return random_energies

