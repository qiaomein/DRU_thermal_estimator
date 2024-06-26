#%%

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.signal as sig
import control as ct
import csv
from datetime import datetime
import re
from scipy.signal import savgol_filter

from parameters import *



test_data_filename = "gel_test_data.csv"

#%%

def f2k(f):
    return (f-32)/1.8 + 273.15


def k2f(k):
    return (k-273.15) * 1.8 + 32


class TimeStamp(object):
    def __init__(self):
        pass



class MysteriousFluid(object):
    def __init__(self,name) -> None:
        self.name = name
        self.timestamps, self.timespan, self.heater_on, self.temp_cell, self.temp_fluid, self.dial_readings = [],[],[],[],[],[]
        self.time_dict = dict()

    def addRawData(self,row):
        timestamp,heater,celltemp,fluidtemp,dialreading = [row[0],row[4],row[8],row[9],row[10]]

        if len(self.timestamps) == 0:
            print("FLUID DATA REFERENCED AT ", timestamp)
            tnum = 0
        else:
            tnum = (self.processTimestamp(timestamp))

        self.timespan.append(tnum)


        self.timestamps.append(timestamp)
        self.time_dict[int(tnum)] = timestamp
        
        self.heater_on.append(heater)
        self.temp_cell.append(f2k(float(celltemp)))
        self.temp_fluid.append(f2k(float(fluidtemp)))
        self.dial_readings.append(float(dialreading))

    def finalizeRawData(self):
        self.timespan = np.array(self.timespan)
        self.temp_cell = np.array(self.temp_cell)
        self.temp_fluid = np.array(self.temp_fluid)
        self.heater_on = np.array(self.heater_on)
        self.dial_readings = np.array(self.dial_readings)

    def processTimestamp(self,timestamp):
        assert self.timestamps
        d = re.split('\ |-|:',timestamp)
        pd = re.split('\ |-|:',self.timestamps[-1])
        currd = datetime(int('20'+d[2]),6,int(d[0]), int(d[3]),int(d[4]),int(d[5]),int(d[6]))
        prevd = datetime(int('20'+pd[2]),6,int(pd[0]), int(pd[3]),int(pd[4]),int(pd[5]),int(pd[6]))
        elapsed = currd - prevd
        

        return elapsed.total_seconds() + self.timespan[-1]


def simple_model(t, Tfi, Tcell, alpha, hgap = heater_gap):
    
    return np.exp(-alpha*cell_area*t/(volume_f*hgap)) * (Tfi - Tcell) + Tcell

#%%
if __name__ == "__main__": # KEEPING ALL TEMPS IN KELVIN

    # read in test data
    my_fluid = MysteriousFluid("testfluid_1")


    with open(test_data_filename, newline ='') as csvfile:
        creader = csv.reader(csvfile)
        next(creader) # skip header
        c = 0
        for row in creader:
            c += 1
            my_fluid.addRawData(row)
            if c > 100000:
                break

    my_fluid.finalizeRawData()

    print("volume of fluid: ", volume_f*1e6, ' mL')


    t1 = 9100
    t2 = t1+800

    t1 = 65300
    t2 = t1 + 3600*2

    # do simple model experiment
    for alpha in [1e-9, 1e-7,3e-7, 1e-6, 1e-5]:

        Tfi = 318.8 #initial fluid temperature in K
        Tcell = 322.8 # cell temp in K

        Tfi = 319.5
        Tcell = 323

        tspan = np.linspace(0,t2-t1,1000)
        ftemps = simple_model(tspan,Tfi,Tcell, alpha, hgap=.018)

        

        
        indices = ( my_fluid.timespan < t2 ) & (my_fluid.timespan > t1)
        my_fluid_timespan = my_fluid.timespan[indices]
        my_fluid_timespan -= my_fluid_timespan[0]





        my_fluid_temps = savgol_filter(my_fluid.temp_fluid[indices],10,5) # smooth out some of the data
        my_fluid_celltemps = my_fluid.temp_cell[indices]
        my_fluid_power = my_fluid.heater_on[indices]


        plt.figure(figsize=(15,10))
        plt.subplot(411)
        # plt.hlines(f2k(120),min(tspan),max(tspan))
        plt.plot(tspan, ftemps)
        plt.plot(my_fluid_timespan,my_fluid_temps)
        plt.tick_params('x', labelbottom=False)
        plt.ylabel("temperature [K]")


        # look in log space

        Tratio = (ftemps - Tcell) / (Tfi - Tcell)
        my_fluid_Tratio = (my_fluid_temps - Tcell) / (Tfi - Tcell)

        plt.subplot(412)
        plt.semilogy(tspan,Tratio)
        plt.semilogy(my_fluid_timespan,my_fluid_Tratio)
        plt.tick_params('x', labelbottom=False)
        plt.ylabel("(Tf-Tc) / (Tfi - Tc)[K/K]")

        # plot cell temperature
        plt.subplot(413)
        plt.plot(tspan,Tcell*np.ones((len(tspan),1)))
        plt.plot(my_fluid_timespan,my_fluid_celltemps)
        
        plt.xlabel("time [s]")
        #plt.tick_params('x', labelbottom=False)
        plt.ylabel("Cell Temperature [K]")

        
        """# plot heater
        plt.subplot(414)
        plt.plot(my_fluid_timespan,my_fluid_power)
        plt.xlabel("time [s]")
        plt.ylabel("Heater on [0/1]")"""


        plt.suptitle(f"Specific heat of mud [J/kg K]: {0.25/(alpha*density_f)}")





# %%
