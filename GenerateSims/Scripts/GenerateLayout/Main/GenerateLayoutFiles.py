#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 14:12:38 2024

@author: chiche
"""

from IceParameters import *
from GenerateStarshapes import generate_arm
from FunctionsGenerateLayout import generate_input_file


# =============================================================================
#                       Antenna positions generation
# =============================================================================

# Starshape layout
x_pos, y_pos, z_pos = \
generate_arm(AntArm, step, angles, IceLevel, AntennaDepth)

x_pos, y_pos, z_pos = x_pos, y_pos, z_pos

# =============================================================================
#                       Input files building
# =============================================================================

# Maximal number antennas per group
NantMax = 5
SimDir = "TestPNFS"
        

generate_input_file(x_pos, y_pos, z_pos, NantMax, SimDir)

