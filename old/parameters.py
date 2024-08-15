"""
EVERYTHING IS IN SI UNITS
"""

import numpy as np

# heater params
power = 250 # in Watts
d_heat = 0.25*0.0254 # m
l_heat = 2*0.0254 # m
a_heat = np.pi*d_heat**2/4*l_heat # m2
n_heaters = 4

# cell geometry params (SS 316)
cap_ss = 502 #in J/kgK
cell_area = 0.0226 # from solidworks; area of cell touching fluid
heater_gap = 0.0057 
k_ss = 18.9 # in W/mK
density_ss = 8000 # in kg/m3
volume_ss = 3871.878*1e-6 * .13
m_ss = density_ss * volume_ss

# fluid properties
cap_f = None
density_f = 12.8 * 119.826427 # converted from ppg -> kg/m3
#volume_f = 0.0729**2/4 * np.pi * 0.1; # this should be constant;
volume_f = 1e-9 * (4.7**2*np.pi * 45.212 + 1/3*24.56*(8.334**2*np.pi + 33.34**2*np.pi + np.sqrt(8.334**2*np.pi*33.34**2*np.pi)) + 33.34**2*np.pi * 76.708 )# all native values are in mm
mass_f = density_f * volume_f
