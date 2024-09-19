#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep 18 19:08:35 2024

@author: chiche
"""


import numpy as np
import matplotlib.pyplot as plt
from scipy.integrate import quad

# Logarithmic scale of energies
E = np.logspace(16, 18, 1000)

# Generating a random energy following a spectrum
def generate_random_energy(CDF, E):
    # Random number between 0 and 1
    u = np.random.rand()
    
    # Finding the index corresping to u in the CDF
    energy_index = np.searchsorted(CDF, u) - 1
    
    # We return the corresponding energy
    return E[energy_index]

# =============================================================================
#                        Test spectrum function
# =============================================================================

def Espectrum(E, gamma):
    return E**(-gamma)
# Spectrum for gamma = 2.7
Espec = Espectrum(E, 2.7)

# Plot of the spectrum
plt.plot(E, Espec)
plt.xscale("log")
plt.yscale("log")
plt.xlabel('Energy (log)')
plt.ylabel('J(E)')
plt.show()

# Total integral for normalization
Ntot = quad(Espectrum, min(E), max(E), args=(2.7,))[0]

# Cummulative Distribution Function
CDF = np.zeros(len(E))

for i in range(1, len(E)):
    # Cummulative integral from Emin to E[i]
    CDF[i] = quad(Espectrum, min(E), E[i], args=(2.7,))[0] / Ntot

# Check for the CDF
plt.plot(E, CDF)
plt.xscale("log")
plt.xlabel('Energy (log)')
plt.ylabel('CDF')
plt.show()

# Generating 10 random energies following the spectrum
random_energies = [generate_random_energy(CDF, E) for _ in range(10000)]

# We check the results
plt.hist(random_energies, bins=np.logspace(16.5, 17.5, 50), \
         density=True, alpha=0.6, label='Energies')
plt.plot(E, Espec / Ntot, label="Spectre théorique $E^{-2.7}$", color="red")
plt.show()

# =============================================================================
#                           Auger Spectrun
# =============================================================================

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

Espec = EspectrumAuger(E)

# Plot of the spectrum
plt.plot(E, Espec)
plt.xscale("log")
plt.yscale("log")
plt.xlabel('Energy (log)')
plt.ylabel('J(E)')
plt.show()


# Total integral for normalization
Ntot = quad(EspectrumAuger, min(E), max(E))[0]

# Cummulative Distribution Function
CDF = np.zeros(len(E))

for i in range(1, len(E)):
    # Cummulative integral from Emin to E[i]
    CDF[i] = quad(EspectrumAuger, min(E), E[i])[0] / Ntot

# Check for the CDF
plt.plot(E, CDF)
plt.xscale("log")
plt.xlabel('Energy (log)')
plt.ylabel('CDF')
plt.show()

# Generating random energies following the spectrum
random_energies = [generate_random_energy(CDF, E) for _ in range(5000)]

Plot = True

if(Plot):
    # We check the results
    plt.hist(random_energies, bins=np.logspace(16, 18, 50), \
             density=True, alpha=0.6, label='Random energies', edgecolor = "black")
    plt.plot(E, Espec / Ntot, label="Auger energy spectrum", color="red")
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("E [eV]")
    plt.ylabel("Counts")
    plt.legend()
    #plt.savefig("E_distrib.pdf", bbox_inches = "tight")
    plt.show()


