
%matplotlib widget

from data_visualizer import MysteriousMud
from parameters import *
import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit

#%%
if __name__ == "__main__": # KEEPING ALL TEMPS IN KELVIN

    # read in test data
    my_fluid = MysteriousMud("water")
    my_fluid.importRawData("rj_data.csv")

    print("volume of fluid: ", volume_f*1e6, ' mL')


    t1 = 0
    t2 = my_fluid.tvector[-1]
    plotting_indices = ( my_fluid.tvector <= t2 ) & (my_fluid.tvector >= t1)
    my_fluid_tvec = my_fluid.tvector[plotting_indices]


    #my_fluid_temps = savgol_filter(my_fluid.temp_fluid[plotting_indices],10,5) # smooth out some of the data
    my_fluid_temps = my_fluid.temp_fluid[plotting_indices]
    #my_fluid_celltemps = savgol_filter(my_fluid.temp_cell[plotting_indices],10,5) # smooth out some of the data
    my_fluid_celltemps = my_fluid.temp_cell[plotting_indices]
    my_fluid_power = my_fluid.heater_on[plotting_indices]
    avg_temp = (my_fluid_celltemps + my_fluid_temps)/2
    davg_temp = avg_temp-Tamb
    total_C = mass_ss*cap_ss + mass_f*4182
    T0 = avg_temp[0]

    decay_model = lambda t,R: (T0-Tamb)*np.exp(-t/(R*total_C))
    
    Rj_estimate, cov = curve_fit(decay_model,my_fluid_tvec,davg_temp)
    Rj_list = np.array([Rj_estimate*.8, Rj_estimate, Rj_estimate*1.2])
    decay_plot = decay_model(my_fluid_tvec,Rj_list)
    print(Rj_estimate, cov)


    ### plotting below

    plt.figure()
 
    # plt.hlines(f2k(120),min(tspan),max(tspan))
    plt.plot(my_fluid_tvec, my_fluid_temps)
    plt.plot(my_fluid_tvec,my_fluid_celltemps)

    plt.ylabel("temperature [K]")


    plt.figure()
 
    # plt.hlines(f2k(120),min(tspan),max(tspan))
    plt.plot(my_fluid_tvec, my_fluid_temps-my_fluid_celltemps)

    plt.ylabel("delta temperature [K]")


    plt.figure()
    plt.plot(my_fluid_tvec,davg_temp)
    for i in range(3):
        plt.plot(my_fluid_tvec,decay_plot[i,:],label = f"Rj = {np.round(Rj_list[i],3)})")
    #plt.title(np.round(Rj_list,4))
    plt.legend()




