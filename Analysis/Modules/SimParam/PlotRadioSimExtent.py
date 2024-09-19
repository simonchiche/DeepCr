#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Sep  4 22:31:01 2024

@author: chiche
"""

import numpy as np
import matplotlib.pyplot as plt

def PlotRadioSimExtent(Depths, radioextent, simextent):
    
    # radio extent and simulation extent vs depth
    plt.plot(3216-np.array(Depths), radioextent, label = "radio")
    plt.plot(3216-np.array(Depths), simextent, label = "sim")
    plt.xlabel("Depth [m]")
    plt.ylabel("Extent [m]")
    #plt.title("E =$%.2f\,$EeV, $\\theta=%.1f^{\circ}$" %(energy, theta), size =14)
    plt.legend()
    #plt.savefig(OutputPath + "LayoutExtent_E%.2f_th%.1f.pdf" \
     #            %(energy, theta), bbox_inches = "tight")
    plt.show()



def PlotFillingFactor(Depths, radioextent, simextent):

    # filling factor
    plt.scatter(3216-np.array(Depths), radioextent/simextent)
    plt.xlabel("Depth [m]")
    plt.ylabel("Filling factor [%]")
    #plt.title("E =$%.2f\,$EeV, $\\theta=%.1f^{\circ}$" %(energy, theta), size =14)
    #plt.savefig(OutputPath + "FillingFactor_E%.2f_th%.1f.pdf" \
    #             %(energy, theta), bbox_inches = "tight")
    plt.show()