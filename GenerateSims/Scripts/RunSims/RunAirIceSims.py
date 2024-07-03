#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 16 05:17:26 2024

@author: chiche
"""
# ssh -X -o ServerAliveInterval=100 schiche@mshort.iihe.ac.be


import numpy as np
import subprocess
import glob
import os

HostDir = "/user/schiche/AirIceSims/InputFiles/Spole1Exp"
Simulations =glob.glob(HostDir)
step = 2


for i in range(len(Simulations)):
    
    filename = Simulations[i]

    if(step==2):
        cmd2_Coreas =  'condor_submit /user/schiche/AirIceSims/MakeSimFiles/job_files\
        /2_Coreas/run_corsika.submit \
        INPUT_FILE="/user/schiche/AirIceSims/InputFiles/Stshp_Proton_0.01_0.0_0.0_1\
        /'+ filename +'.inp" RUN_NR="000001" ATMOS_FILE="/user/schiche\
        /AirIceSims/InputFiles/TestSim1/Atmosphere.dat" REAS_FILE="/user/schiche\
        /AirIceSims/InputFiles/TestSim1/SIM.reas" LIST_FILE="/user/schiche/AirIceSims\
        /InputFiles/TestSim1/SIM.list" OUTPUT_DIR="/user/schiche/AirIceSims/OutputFiles\
        /CoreasOutput/'+ filename +'" LOG_NAME="="'
        p = subprocess.Popen(cmd2_Coreas, cwd=os.getcwd(), shell=True)
        stdout, stderr = p.communicate()
        
    
    if(step==3):
        cmd3_GeantInput = 'condor_submit make_job_files.submit INPUT_FILE="\
        /user/schiche/AirIceSims/OutputFiles/CoreasOutput/Stshp_Proton_0.01_0.0\
        _0.0_1/DAT000001" OUTPUT_DIR="/user/schiche/AirIceSims/OutputFiles/\
        GeantInputFiles/Stshp_Proton_0.01_0.0_0.0_1" ZENITH="0.0" AZIMUTH="0.0" \
        CORSIKA_LOG_FILE="/user/schiche/AirIceSims/OutputFiles/CoreasOutput/\
        Stshp_Proton_0.01_0.0_0.0_1/RUN000001.log" LOG_NAME="Stshp_Proton_\
        0.01_0.0_0.0_1"'
        p = subprocess.Popen(cmd3_GeantInput, cwd=os.getcwd(), shell=True)
        stdout, stderr = p.communicate()
        
        
    if(step ==4):
        cmd4_Dagfile = 'python /user/schiche/AirIceSims/MakeSimFiles\
        /python_scripts/4_GeantSims/make_dagfile.py /user/schiche/AirIceSims\
        /OutputFiles/GeantInputFiles/Stshp_Proton_0.01_0.0_0.0_1 /user/schiche\
        /AirIceSims/OutputFiles/GeantOutput/Stshp_Proton_0.01_0.0_0.0_1 \
        Stshp_Proton_0.01_0.0_0.0_1 /user/schiche/AirIceSims/OutputFiles\
        /DagOutput/Stshp_Proton_0.01_0.0_0.0_1 0.0 0.0 /user/schiche/AirIceSims\
        /InputFiles/Stshp_Proton_0.01_0.0_0.0_1/SIM.reas /user/schiche/AirIceSims\
        /InputFiles/Stshp_Proton_0.01_0.0_0.0_1/SIM.list /user/schiche/AirIceSims\
        /InputFiles/Stshp_Proton_0.01_0.0_0.0_1/Atmosphere.dat'
        p = subprocess.Popen(cmd4_Dagfile, cwd=os.getcwd(), shell=True)
        stdout, stderr = p.communicate()
        
        
    if(step == 4.5):
        cmd4bis_GeantSims = 'condor_submit_dag /user/schiche/AirIceSims/\
        OutputFiles/DagOutput/Stshp_Proton_0.01_0.0_0.0_1/run_ice_shelf_project.dag'
        p = subprocess.Popen(cmd4bis_GeantSims, cwd=os.getcwd(), shell=True)
        stdout, stderr = p.communicate()
        
        
    if(step ==5):
        cmd5_CheckSims = 'python3 /user/schiche/AirIceSims/OutputFiles\
        /GeantInputFiles/Stshp_Proton_0.01_0.0_0.0_1/check_for_killed_jobs.py \
        /user/schiche/AirIceSims/OutputFiles/DagOutput/Stshp_Proton_0.01_0.0_0.\
        0_1/run_ice_shelf_project.dag'
        p = subprocess.Popen(cmd5_CheckSims, cwd=os.getcwd(), shell=True)
        stdout, stderr = p.communicate()
        
    
    if(step == 6): 
        cmd6_CombineOutputs = 'condor_submit /user/schiche/AirIceSims\
        /MakeSimFiles/job_files/6_CombineGeantOutputs/run_combine_root_\
        files.submit INPUT_DIR="/user/schiche/AirIceSims/OutputFiles\
        /GeantOutput/Stshp_Proton_0.01_0.0_0.0_1" LOG_NAME="Stshp_Proton_0.01_0.0_0.0_1"'
        p = subprocess.Popen(cmd6_CombineOutputs, cwd=os.getcwd(), shell=True)
        stdout, stderr = p.communicate()
        
        
        cmd6bis = 'condor_submit /user/schiche/AirIceSims/MakeSimFiles\
        /job_files/6_CombineGeantOutputs/run_combine_antenna_traces.submit \
        INPUT_DIR="/user/schiche/AirIceSims/OutputFiles/GeantOutput\
        /Stshp_Proton_0.01_0.0_0.0_1" LOG_NAME="Stshp_Proton_0.01_0.0_0.0_1"'
        p = subprocess.Popen(cmd6bis, cwd=os.getcwd(), shell=True)
        stdout, stderr = p.communicate()

    

        
        cmd2_Coreas = 'condor_submit /user/schiche/AirIceSims/MakeSimFiles/job_files\
        /2_Coreas/run_corsika.submit \
        INPUT_FILE='+ SimfromDir + '\
        /'+ filename +'.inp" RUN_NR="000001" ATMOS_FILE="'+ SimfromDir + '/Atmosphere.dat" REAS_FILE='+ SimfromDir + '/SIM.reas" LIST_FILE="/user/schiche/AirIceSims\
        #/InputFiles' + SimfromDir + '/SIM.list" OUTPUT_DIR="/user/schiche/AirIceSims/OutputFiles\
        #/CoreasOutput/'+ filename +'" LOG_NAME="="'