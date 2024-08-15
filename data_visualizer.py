#%%

#%matplotlib widget

import numpy as np
import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
from scipy.optimize import curve_fit
import scipy.signal as sig
import csv
from datetime import datetime
import re


from parameters import *



test_data_filename = "gel_test_data.csv"
sys_data_filename = "systrend.csv"

#%%

def f2k(f):
    return (f-32)/1.8 + 273.15


def k2f(k):
    return (k-273.15) * 1.8 + 32


class Timestamp(object):
    def __init__(self,raw_timestamps): # raw_timestamps is a list of the raw timestamps
        self.raw_timestamps = raw_timestamps
        self._month_dict = dict()
        months = "Jan,Feb,Mar,Apr,May,Jun,Jul,Aug,Sep,Oct,Nov,Dec".split(',')
        for i,m in enumerate(months):
            self._month_dict[m] = i+1

        self.timespan = self.convertRawTimestamps(raw_timestamps)


    def convertRawTimestamps(self, all_timestamps):
        timespan = [0]
        N = len(all_timestamps)
        for i in range(1,N):
            d = re.split('\ |-|:',all_timestamps[i]) # current date
            pd = re.split('\ |-|:',all_timestamps[i-1]) #previous date
            #print(d,pd)
            n_month = self._month_dict[d[1]]
            pn_month = self._month_dict[pd[1]]

            currd = datetime(int('20'+d[2]), n_month, int(d[0]), int(d[3]),int(d[4]),int(d[5]),int(d[6]))
            prevd = datetime(int('20'+pd[2]), pn_month, int(pd[0]), int(pd[3]),int(pd[4]),int(pd[5]),int(pd[6]))
            
            #print(currd, prevd)
            
            elapsed = currd - prevd
            elapsed_t = elapsed.total_seconds()
            if elapsed_t > 10:
                print(f"uh oh... something went horribly wrong! at dates {currd, prevd}")
            timespan.append(elapsed_t + timespan[-1])
        
        return timespan

    #TODO: add a method to handle multiple timestamps with close but not necessarily same values


class MysteriousMud(object):
    def __init__(self,name) -> None:
        self.name = name
        self.tvector, self.heater_on, self.temp_cell, self.temp_fluid, self.dial_readings = [],[],[],[],[]
        self.timestamp = None
        self.MAX_ENTRIES = 1e6
        self.filename = None
        self.ref_timestamp = None
        

    def importRawData(self,filename, heaterHeader = True): #heaterHeader is false if using old data (no GELHtrOut(%) header field); if there is that field, this flag is true
        
        """
        the importRawData method takes in the filename of the csv file with raw data and imports it into the object
        """
        
        full_dir = f"data\\{filename}"
        timestamps  = []
        self.filename = filename
        try:

            with open(full_dir, newline ='') as csvfile:
                creader = csv.reader(csvfile)
                next(creader) # skip header
                for c,row in enumerate(creader): # loop through each row in csvfile
                    if c > self.MAX_ENTRIES:
                        print("TOO MANY ENTIRES!")
                        break
                    
                    if heaterHeader:
                        timestamp,heater,celltemp,fluidtemp,dialreading = [row[0],row[7],row[8],row[9],row[10]]
                    else:
                        timestamp,heater,celltemp,fluidtemp,dialreading = [row[0],row[4],row[8],row[9],row[10]]
                    
                    if not timestamps: # if timestamps is empty
                        print(f"FLUID DATA FOR {self.name} REFERENCED AT ", timestamp)
                        self.ref_timestamp = timestamp
                    timestamps.append(timestamp)
                    
                    self.heater_on.append(float(heater))
                    self.temp_cell.append(f2k(float(celltemp)))
                    self.temp_fluid.append(f2k(float(fluidtemp)))
                    self.dial_readings.append(float(dialreading))

        except FileNotFoundError:
            print(f"{filename} file not found!")
            return None

        self.timestamp = Timestamp(timestamps)

        self.tvector = np.array(self.timestamp.timespan)
        self.temp_cell = np.array(self.temp_cell)
        self.temp_fluid = np.array(self.temp_fluid)
        self.heater_on = np.array(self.heater_on)
        if heaterHeader:
            self.heater_dutycycle = self.heater_on
        else:
            self.heater_dutycycle = None
        self.dial_readings = np.array(self.dial_readings)


#%%
if __name__ == "__main__": # KEEPING ALL TEMPS IN KELVIN

    # this main script is to test the MysteriousMud class and Timestamp classes on example data

    # read in test data
    my_fluid = MysteriousMud("testfluid_1")
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
    Rj_estimate,_ = curve_fit(decay_model,my_fluid_tvec,davg_temp)
    decay_plot = decay_model(my_fluid_tvec,Rj_estimate)


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
    plt.semilogy(my_fluid_tvec,davg_temp)
    plt.semilogy(my_fluid_tvec,decay_plot)


