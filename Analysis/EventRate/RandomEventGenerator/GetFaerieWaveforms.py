#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 14:38:32 2024

@author: chiche
"""

import sys
#sys.path.append('/Users/chiche/Library/Python/3.9/lib/python/site-packages')
sys.path.append("/Users/chiche/Desktop/NuRadioMC")
sys.path.append("/Users/chiche/Desktop/NuRadioMC/NuRadioReco")

import NuRadioMC
import NuRadioReco.modules.io.coreas.readFEARIEShower

from NuRadioReco.modules import efieldToVoltageConverter

from NuRadioReco.detector import generic_detector as detector

from NuRadioReco.utilities import units

import argparse
from matplotlib import pyplot as plt
import sys
import numpy as np

def GetFaerieVoltage(inputfilename, detectordescription):

    det = detector.GenericDetector(detectordescription)
    
    efield_converter = efieldToVoltageConverter.efieldToVoltageConverter()
    efield_converter.begin()

    readFEARIEShower = NuRadioReco.modules.io.coreas.readFEARIEShower.readFEARIEShower()
    readFEARIEShower.begin(inputfilename, det)

    
    for event, det in readFEARIEShower.run(depth=100):
        
        print('Event {} {}'.format(event.get_run_number(), event.get_id()))
        print('Number of stations: {}'.format(len(list(event.get_stations()))))
        
        VoltVpole = dict() 
        tVpole = dict()
        Trigger = []
        
        i =0
        for station in event.get_stations():

            sim_station = station.get_sim_station()

            if len(sim_station.get_electric_fields()) < 2:
                continue

            print(f"Number of electric fields: {len(sim_station.get_electric_fields())}")
            for efield in sim_station.get_electric_fields():
                print(efield.get_position())

            efield_converter.run(event, station, det)
            print(f"Number of channels: {len(station.get_channel_ids())}")
            channel = station.get_channel(0)
            
            t = channel.get_times()

            print("time shape:")
            print(t.shape)
            V = channel.get_trace() / units.mV
            VoltVpole[i] = V
            tVpole[i] = t
            Trigger.append(event.has_triggered())
            print("Channel shape:")
            print(np.shape(V))
            i = i +1

    return VoltVpole, tVpole, Trigger

inputfilename = ["/Users/chiche/Desktop/REGevents/HDF5files/Rectangle_Proton_0.1_10_0_1_0.hdf5"]
detectordescription = "/Users/chiche/Desktop/DeepCrSearch/NuRadioMC/NuRadioReco/detector/RNO_G/RNO_single_station.json"


VoltVpole, tVpole, Trigger = GetFaerieVoltage(inputfilename, detectordescription)


print(len(VoltVpole))