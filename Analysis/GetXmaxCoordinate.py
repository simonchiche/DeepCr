#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 19 17:01:37 2024

@author: chiche
"""


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Tue Apr 30 16:47:26 2024

@author: chiche
"""
import numpy as np

def showerdirection(azimuth, zenith_cr):
    
    # zenith in degrees, cosmic-ray convention
    # azimuth in degrees (0: shower towards north, 90: towards east)
    
    zenith = 180 - zenith_cr

    zenith = zenith*np.pi/180.0
    azimuth = azimuth*np.pi/180.0
    
    uv = np.array([np.sin(zenith)*np.cos(azimuth), \
                   np.sin(zenith)*np.sin(azimuth), np.cos(zenith)])
    
    
    return uv

#print(showerdirection(0, 150)*100)

#print(showerdirection(270, 117))

# check zenith convention

def _getAirDensity(_height, model):

    '''Returns the air density at a specific height, using either an 
    isothermal model or the Linsley atmoshperic model as in ZHAireS

    Parameters:
    ---------
        h: float
            height in meters

    Returns:
    -------
        rho: float
            air density in g/cm3
    '''

    if model == "isothermal":
            #Using isothermal Model
            rho_0 = 1.225*0.001    #kg/m^3
            M = 0.028966    #kg/mol
            g = 9.81        #m.s^-2
            T = 288.        #
            R = 8.32        #J/K/mol , J=kg m2/s2
            rho = rho_0*np.exp(-g*M*_height/(R*T))  # kg/m3

    elif model == "linsley":
        #Using Linsey's Model
        bl = np.array([1222., 1144., 1305.5948, 540.1778,1])*10  # g/cm2  ==> kg/cm3
        cl = np.array([9941.8638, 8781.5355, 6361.4304, 7721.7016, 1e7])  #m
        hl = np.array([4,10,40,100,113])*1e3  #m
        if (_height>=hl[-1]):  # no more air
            rho = 0
        else:
            hlinf = np.array([0, 4,10,40,100])*1e3  #m
            ind = np.logical_and([_height>=hlinf],[_height<hl])[0]
            rho = bl[ind]/cl[ind]*np.exp(-_height/cl[ind])
            rho = rho[0]*0.001
    else:
        print("#### Error in GetDensity: model can only be isothermal or linsley.")
        return 0

    return rho

_getAirDensity(1000, "linsley")*1e3 # in k.m-3 with the 1e3

def _get_CRzenith(zen, GdAlt, injh):
    ''' Corrects the zenith angle for CR respecting Earth curvature, zenith seen by observer
        ---fix for CR (zenith computed @ shower core position
    
    Arguments:
    ----------
    zen: float
        GRAND zenith in deg
    injh: float
        injection height wrt to sealevel in m
    GdAlt: float
        ground altitude of array/observer in m (should be substituted)
    
    Returns:
    --------
    zen_inj: float
        GRAND zenith computed at shower core position in deg
        
    Note: To be included in other functions   
    '''
    Re= 6370949 # m, Earth radius
    
    a = np.sqrt((Re + injh)**2. - (Re+GdAlt)**2 *np.sin(np.pi-np.deg2rad(zen))**2) - (Re+GdAlt)*np.cos(np.pi-np.deg2rad(zen))
    zen_inj = np.rad2deg(np.pi-np.arccos((a**2 +(Re+injh)**2 -Re**2)/(2*a*(Re+injh))))
    
    return zen_inj    

        
def _dist_decay_Xmax(Xmax_primary, gdalt, zenith_cr, injh2): 
    ''' Calculate the height of Xmax and the distance injection point to Xmax along the shower axis
    
    Arguments:
    ----------
    zen: float
        CR zenith in deg, for CR shower use _get_CRzenith()
    injh2: float
        injectionheight above sealevel in m
    Xmax_primary: float
        Xmax in g/cm2 
        
    Returns:
    --------
    h: float
        vertical Xmax_height in m
    ai: float
        Xmax_distance injection to Xmax along shower axis in m
    '''
     
    zenith_grand = 180 - zenith_cr 
    if(zenith_cr > 40):
        zenith_grand = _get_CRzenith(zenith_grand, gdalt, injh2)

    zen2 = np.deg2rad(zenith_grand)
    
    hD=injh2
    step=10 #m
    if hD>10000:
        step=10 #m
    Xmax_primary= Xmax_primary#* 10. # g/cm2 to kg/m2: 1g/cm2 = 10kg/m2
    gamma=np.pi-zen2 # counterpart of where it goes to
    Re= 6370949 # m, Earth radius
    X=0.
    i=0.
    h=hD
    ai=0
    while X< Xmax_primary:
        i=i+1
        ai=i*step #100. #m
        hi= -Re+np.sqrt(Re**2. + ai**2. + hD**2. + 2.*Re*hD - 2*ai*np.cos(gamma) *(Re+hD))## cos(gamma)= + to - at 90dg
        deltah= abs(h-hi) #(h_i-1 - hi)= delta h
        h=hi # new height
        rho = _getAirDensity(hi, "linsley") ## Nicolas you can put your air density model here
        
        X=X+ rho * step*100. #(deltah*100) *abs(1./np.cos(np.pi-zen2)) # Xmax in g/cm2, slanted = Xmax, vertical/ cos(theta); density in g/cm3, h: m->100cm, np.pi-zen2 since it is defined as where the showers comes from, abs(cosine) so correct for minus values
       
    return h, ai # Xmax_height in m, Xmax_distance in m    

def getGroundXmaxDistance(Xmax_primary, zenith_cr, injh2, GroundAltitude):
    

    # CR zenith here? check with RM
    XmaxHeight, DistDecayXmax = _dist_decay_Xmax(Xmax_primary, GroundAltitude, zenith_cr, injh2) 

    print(XmaxHeight, "height")
    Rearth = 6370949 
    
    # zenith in cosmic ray convention here
    
    zenith_cr = zenith_cr*np.pi/180
    
    dist = np.sqrt((Rearth+ XmaxHeight)**2 - ((Rearth + GroundAltitude)*np.sin(zenith_cr))**2) \
    - (Rearth + GroundAltitude)*np.cos(zenith_cr)
    
    if(XmaxHeight< GroundAltitude): dist =0
            
    return dist   

def getXmaxPosition(XmaxPrimary, azimuth, zenith_cr, glevel,  injh=1e5):
    
    # zenith in CR convention
    # azimuth in ZHS Coreas conv
    uv = showerdirection(azimuth, zenith_cr)
    
    showerDistance = getGroundXmaxDistance(XmaxPrimary, zenith_cr, injh, glevel)
    XmaxPosition = -uv*showerDistance 
    XmaxPosition[2] = XmaxPosition[2] + glevel  
    
    return XmaxPosition

getXmaxPosition(765.70844, 0, 70, 3216)


    