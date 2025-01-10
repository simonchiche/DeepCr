#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 19 14:38:32 2024

@author: chiche
"""

import sys
import logging
import matplotlib.pyplot as plt
import NuRadioReco.modules.channelGenericNoiseAdder
import NuRadioReco.modules.RNO_G.hardwareResponseIncorporator
import NuRadioReco.modules.trigger.highLowThreshold as highLowThreshold
import NuRadioReco.modules.trigger.powerIntegration as powerIntegration
import NuRadioReco.modules.trigger.simpleThreshold as simpleThreshold
import NuRadioReco.modules.channelBandPassFilter
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
from icecream import ic

#logger = logging.getLogger()
#logger.setLevel(logging.CRITICAL)

channelGenericNoiseAdder = NuRadioReco.modules.channelGenericNoiseAdder.channelGenericNoiseAdder()
channelGenericNoiseAdder.begin()
hardwareResponseIncorporator = NuRadioReco.modules.RNO_G.hardwareResponseIncorporator.hardwareResponseIncorporator()

channelBandPassFilter = NuRadioReco.modules.channelBandPassFilter.channelBandPassFilter()
channelBandPassFilter.begin()   

def GetFaerieVoltage(inputfilename, detectordescription):

    """
    Take an input hdf5 file and a detector description and run the FAERIE reader to return waveforms
    inputfilename: hdf5 file generated with the coreas converter (geant observers option enabled)
    detectordescription: json file describing each observer
    """

    det = detector.GenericDetector(detectordescription)
    ic('Initalizing Efield to Volatge Conversion')
    efield_converter = efieldToVoltageConverter.efieldToVoltageConverter()
    efield_converter.begin()
    ic('Initalizing FAERIE reader')
    readFEARIEShower = NuRadioReco.modules.io.coreas.readFEARIEShower.readFEARIEShower()
    readFEARIEShower.begin(inputfilename, det)
    ic(inputfilename, det)

    ic('Loop over the events, running FAERIE reader')
    for event, det in readFEARIEShower.run(depth=100):
        ic('we are here')
        print('Event {} {}'.format(event.get_run_number(), event.get_id()))
        print('Number of stations: {}'.format(len(list(event.get_stations()))))
        
        #### INITIALISATION #####
        VoltVpole = dict() 
        tVpole = dict()
        freqChannel = dict()
        SpectrumChannel = dict()
        Antpos = []
        Trigger =[]
        i =0

        #input("Looping over the stations")
        ic('Loop over the stations')
        for station in event.get_stations():
            
            if(i%100==0):ic('station', i)
            
            ### Getting the stations
            sim_station = station.get_sim_station()
            
            ### Cross checks
            #if len(sim_station.get_electric_fields()) < 2:
            #    continue

            #print(f"Number of electric fields: {len(sim_station.get_electric_fields())}")
            #for efield in sim_station.get_electric_fields():
            #    print(efield.get_position())
            
            try:
                #ic('Efield to voltage conversion')
                efield_converter.run(event, station, det)
            except TypeError:
                print("A TypeError occurred, handling it!")
                freqChannel[i] = np.nan
                SpectrumChannel[i] = np.nan
                VoltVpole[i] = np.nan
                tVpole[i] = np.nan
                #Trigger.append(np.nan)
                #i = i +1
                continue

            for efield in sim_station.get_electric_fields():
                Antpos.append(efield.get_position())
                break

            ### RF chain (hardware response)
            #ic('Hardware response')

            RFchain = True
            if(RFchain):
                hardwareResponseIncorporator.run(event, station, det, sim_to_data=True)
            
            ##channel = station.get_channel(0)
            ##VoltageTrace = channel.get_trace() / units.mV
            ##trace = channel.get_trace()
            ##zero_trace = np.zeros_like(trace)
            ##channel.set_trace(zero_trace, 0.5)

            #ic('Noise generation')
            Noise = True
            if(Noise):
                NoiseLevel = 4 * units.mV
                channelGenericNoiseAdder.run(event, station, det,
                                    amplitude=NoiseLevel,#cfg['Vrms_thermal_noise'],
                                    min_freq=50 * units.MHz, #cfg['T_noise_min_freq'],  # check frequency band
                                    max_freq=2000 * units.MHz, #cfg['T_noise_max_freq'], 
                                    type='rayleigh', bandwidth=None)
                
            channel = station.get_channel(0)
            VoltageTrace = channel.get_trace() #* units.mV
            #ic(max(abs(VoltageTrace)))
         
            ### Getting chanel info
            #print(f"Number of channels: {len(station.get_channel_ids())}")
            channel = station.get_channel(0)
            TimeTrace = channel.get_times()
            freqChannel[i] = channel.get_frequencies()
            SpectrumChannel[i] = np.abs(channel.get_frequency_spectrum())

            ### Trigger methods
            threshold = 5*NoiseLevel
            #print(threshold)
            dt = 5 * units.ns
            #triggered_samples = highLowThreshold.get_high_low_triggers(VoltageTrace, threshold, -threshold, dt)
            #triggered_samples = highLowThreshold.get_high_low_triggers(VoltageTrace, threshold, -threshold, dt)
            triggered_samples =  simpleThreshold.get_threshold_triggers(VoltageTrace, threshold)

            if True in triggered_samples:
                Trigger.append(True)
            else:
                Trigger.append(False)
            VoltVpole[i] = VoltageTrace
            tVpole[i] = TimeTrace
            
            i = i +1
            
            Plot = False
            if(Plot):
                fig, axs = plt.subplots(1, 2)
                axs[0].plot(channel.get_times(), channel.get_trace() / units.mV)
        
                axs[0].set_xlabel('time / ns')
                axs[0].set_ylabel('voltage / mV')
        
                axs[1].plot(channel.get_frequencies(), np.abs(channel.get_frequency_spectrum()))
                axs[1].set_xlabel('frequencies / GHz')
                axs[1].set_xlim(None, 1.2)
        
                plt.show()

    return VoltVpole, tVpole, Trigger, threshold, freqChannel, SpectrumChannel, np.array(Antpos)
""""
### Running the data processor
inputfilename = ["/Users/chiche/Desktop/REGevents/HDF5files/Rectangle_Proton_0.1_10_0_1_0.hdf5"]
#inputfilename = [ "/Users/chiche/Desktop/Test/Rectangle_Proton_0.001_0_0_1/Rectangle_Proton_0.001_0_0_1.hdf5"]
detectordescription = "/Users/chiche/Desktop/DeepCrSearch/NuRadioMC/NuRadioReco/detector/RNO_G/RNO_single_channel.json"

VoltVpole, tVpole, Trigger, threshold,  freqChannel, SpectrumChannel, Antpos\
    = GetFaerieVoltage(inputfilename, detectordescription)

### Trigger results
#print(Trigger)
print(sum(Trigger), len(Trigger), sum(Trigger)/len(Trigger))
print(len(VoltVpole))
print(max(abs(VoltVpole[0])))

plt.scatter(Antpos[:,0], Antpos[:,1], label = "All events")
plt.scatter(Antpos[:,0][Trigger], Antpos[:,1][Trigger], color ="red", label = "Triggered events")
plt.xlabel("x [m]")
plt.ylabel("y [m]")
plt.legend()
#plt.savefig("/Users/chiche/Desktop/Trigger.pdf", bbox_inches ="tight")
plt.show()
#sys.exit()
### Voltage Map
Vmax = np.array([np.max(np.abs(arr)) for arr in VoltVpole.values()])

plt.scatter(Antpos[:,0], Antpos[:,1], c=  Vmax, cmap = "jet", label = "All events")
plt.show()

"""
"""
### Analysing the outputs
AmplitudeThreshold = 0
k =0
ilist = np.array([203, 237, 264, 265])
for i in range(len(tVpole)):
   print(max(abs(VoltVpole[i])))
   if(max(abs(VoltVpole[i])>=AmplitudeThreshold) & (k<40)):
   #if(i in ilist):
        print(i)
        ic(np.std(VoltVpole[i]))
        plt.plot(tVpole[i], VoltVpole[i])
        plt.title("%.d" %i)
        plt.xlabel('time / ns')
        plt.ylabel('voltage / mV')
        #plt.savefig("/Users/chiche/Desktop/DeepCrPlots_13_11_24/VoltageNoisy/VoltageTrace_%.d.pdf" %i, bbox_inches = "tight")
        plt.show()
        plt.plot(freqChannel[i], SpectrumChannel[i])
        plt.xlabel('frequencies / GHz')
        plt.title("%.d" %i)
        #plt.savefig("/Users/chiche/Desktop/DeepCrPlots_13_11_24/VoltageNoisy/Spectrum_%.d.pdf" %i, bbox_inches = "tight")
        plt.show()
        k = k + 1
"""

