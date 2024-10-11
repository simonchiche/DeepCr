#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 26 17:06:13 2024

@author: chiche
"""

import json
import numpy as np
import matplotlib.pyplot as plt
with open('/Users/chiche/Desktop/RNO_single_station.json', 'r') as file:
    RNO_single_station = json.load(file)
    
Nchannel =24
x_pos = np.zeros(Nchannel)
y_pos = np.zeros(Nchannel)
z_pos = np.zeros(Nchannel)

# Accessing positions for each channel
channels = RNO_single_station['channels']
for channel_id, channel_info in channels.items():
        x_pos[int(channel_id)] = channel_info['ant_position_x']
        y_pos[int(channel_id)] = channel_info['ant_position_y']
        z_pos[int(channel_id)] = channel_info['ant_position_z']  
    

plt.scatter(x_pos, z_pos)
plt.show()