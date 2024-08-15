
%matplotlib widget

from data_visualizer import MysteriousMud, f2k, k2f
from parameters import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

#%%
if __name__ == "__main__": # KEEPING ALL TEMPS IN KELVIN

    # read in test data
    my_fluid = MysteriousMud("water")
    my_fluid.importRawData("step1_10p.csv", heaterHeader=True)

    print("volume of fluid: ", volume_f*1e6, ' mL')


    t1 = my_fluid.tvector[np.where(my_fluid.heater_dutycycle > 0)[0][0]]
    t2 = my_fluid.tvector[-1]
    plotting_indices = ( my_fluid.tvector <= t2 ) & (my_fluid.tvector >= t1)
    my_fluid_tvec = my_fluid.tvector[plotting_indices] - t1


    #my_fluid_temps = savgol_filter(my_fluid.temp_fluid[plotting_indices],10,5) # smooth out some of the data
    my_fluid_temps = my_fluid.temp_fluid[plotting_indices]
    #my_fluid_celltemps = savgol_filter(my_fluid.temp_cell[plotting_indices],10,5) # smooth out some of the data
    my_fluid_celltemps = my_fluid.temp_cell[plotting_indices]
    my_fluid_power = my_fluid.heater_on[plotting_indices]
    avg_temp = my_fluid_temps # 
    
    total_C = mass_ss*cap_ss + mass_f*4182
    T0 = avg_temp[0]
    tempfits = avg_temp-T0

    decay_model = lambda t,T,tau: T*(1-np.exp(-t/(tau)))
    
    transient_indices = np.where(my_fluid_power > 0)[0]
    popt, cov = curve_fit(decay_model,my_fluid_tvec[transient_indices], tempfits[transient_indices], p0 = (40,1000))
    Ta, tau = popt


    ### plotting below

    plt.figure()
 
    # plt.hlines(f2k(120),min(tspan),max(tspan))
    plt.plot(my_fluid_tvec, my_fluid_temps)
    plt.plot(my_fluid_tvec,my_fluid_celltemps)
    plottvec = np.linspace(my_fluid_tvec[0],my_fluid_tvec[-1],1000)
    plt.plot(plottvec, decay_model(plottvec,Ta,tau) + T0)

    plt.ylabel("temperature [K]")
    plt.xlabel("time [s]")

    plt.title(f"{my_fluid.name}: {Ta,tau}")


    plt.figure()
 
    # plt.hlines(f2k(120),min(tspan),max(tspan))
    plt.plot(my_fluid_tvec, my_fluid_power)

    plt.ylabel("heater percentage [%]")




