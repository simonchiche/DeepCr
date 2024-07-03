#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Apr  4 11:27:45 2024

@author: chiche
"""


Energy = 0.5# EeV
TargetE = 1e9 # GeV

def ConvertEnergy(E):
    
    E = E*1e9 # EeV to GeV
    
    Esci = "{:e}".format(E) # scientific notation
    mantissa, exponent = Esci.split('e')
    mantissa = mantissa.rstrip('0').rstrip('.') if '.' in mantissa else mantissa
    exponent = exponent.replace('+', '.') if exponent.startswith('+') else exponent
    Esci = mantissa + 'E' + exponent

    return Esci
    

print(ConvertEnergy(Energy))

