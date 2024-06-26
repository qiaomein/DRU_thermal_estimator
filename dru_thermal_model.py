import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
import scipy.signal as sig
import control as ct
from control.matlab import lsim
import csv
from datetime import datetime
import re
from scipy.signal import savgol_filter
from scipy.integrate import solve_ivp

from parameters import *


class SystemParameters(object):
    def __init__(self) -> None:
        # heater params
        power = 250 # in Watts per heater rod
        d_heat = 0.25*0.0254 # m
        l_heat = 2*0.0254 # m
        a_heat = np.pi*d_heat**2/4*l_heat # m2
        n_heaters = 4

        self.max_heater = n_heaters*power # in W

        # cell geometry params (SS 316)
        self.cap_ss = 502 #in J/kgK
        self.cell_area = 0.0226 # from solidworks; area of cell touching fluid
        self.heater_gap = 0.0057 
        self.k_ss = 18.9 # in W/mK
        self.density_ss = 8000 # in kg/m3
        self.volume_ss = 3871.878*1e-6 * .13
        self.mass_ss = density_ss * volume_ss

        # fluid properties
        self.cap_f = None
        self.density_f = None #12.8 * 119.826427 # converted from ppg -> kg/m3
        #volume_f = 0.0729**2/4 * np.pi * 0.1; # this should be constant;
        self.volume_f = 1e-9 * (4.7**2*np.pi * 45.212 + 1/3*24.56*(8.334**2*np.pi + 33.34**2*np.pi + np.sqrt(8.334**2*np.pi*33.34**2*np.pi)) + 33.34**2*np.pi * 76.708 )# all native values are in mm
        self.mass_f = density_f * volume_f

class ThermalController(object):
    def __init__(self) -> None:
        self.kp = 0
        self.ki = 0
        self.kd = 0
        self.__preverror = 0

        self.reference_temperature = None
    
    def feedback(self,x):
        Tc = x[0]
        Tf = x[1]
        Tref = self.reference_temperature
        
        error = Tref - Tf
        derror = error - self.__preverror
        p = (self.k*(error + self.kd*derror))

        p = self.k*(Tref-Tf) + self.kd*(Tref-Tc) # full state feedback

        self.__preverror = error

        if Tf > Tref:
            p = 0
        
        return p
        

class ThermalModel(object): # one instance of a thermal model (fixed set of parameters)
    def __init__(self, mud) -> None:
        self.mud = mud # this is a MysteriousMud object

        self.parameters = SystemParameters()
        self.controller = ThermalController()
        self.input_power = None
        self.Tambient = mud.temp_cell[0] #300
        self.controller_on = False
        
        # simulation time settings
        self.t_initial = 0
        self.t_final = float(mud.tvector[-1])
        N = int((self.t_final - self.t_initial)//2)
        self.tvec = np.linspace(0,self.t_final,N)

        # data time settings
        

        self.t_stepfinal = self.tvec[np.where(mud.heater_dutycycle > 0)[0][-1]] - self.mud.tvector[np.where(self.mud.heater_dutycycle > 0)[0][0]]



        # parameters
        density_f = 1000 # water
        Rc = heater_gap/(cell_area*k_ss)

        Cf = 4000 * density_f * volume_f # estimate
        Cc = cap_ss * mass_ss
        Rj = 1.86 #0.03/(cell_area*0.3) # thermal resistance of the jacket
        Rf = 1.0428 #.1/(.6*cell_area) #0.02/(cell_area*.4) # estimate equivalent thermal resistance relating to forced convection
        Rfj = 1
        self.thermal_resistances = [Rj,Rfj,Rf]
        self.thermal_capacitances = [Cc, Cf]

    def Tf_final(self, power):
        # final value theorem
        rj,rfj,rf = self.thermal_resistances
        cc,cf = self.thermal_capacitances

        bulk = rfj*rj/((rf+rj)*(rf+rfj))

        return power * bulk

        
    def stepResponse(self, x0 = np.array([0,0]), powerPercentage = 1):
        

        assert powerPercentage > 0

        self.input_power = powerPercentage * self.parameters.max_heater


        # takes in true initial condition x0 = [tc,tf] 
        X0 = x0 - self.Tambient
        sol = solve_ivp(self.__ss_model,(self.t_initial,self.t_final), X0, t_eval=self.tvec)

        yout = sol.y
        self.temp_mud = yout[1,:] + self.Tambient
        self.temp_cell = yout[0,:] + self.Tambient
        print(f"Rj, Rfj, Rf: {self.thermal_resistances}")
        print(f"Cc, Cf: {self.thermal_capacitances}")

        self.final_fluid_temperature = self.Tf_final(self.input_power) + self.Tambient

        return self.temp_cell, self.temp_mud, 
    
    def plot(self, bigfig = True):
        t1 = self.mud.tvector[np.where(self.mud.heater_dutycycle > 0)[0][0]] # start when input is turned on
        t2 = self.mud.tvector[-1]
        self.plotting_indices = ( self.mud.tvector <= t2 ) & (self.mud.tvector >= t1)
        my_fluid_tvec = self.mud.tvector[self.plotting_indices] - t1

        my_fluid_temps = self.mud.temp_fluid[self.plotting_indices]
        my_fluid_celltemps = self.mud.temp_cell[self.plotting_indices]


        if bigfig:
            plt.figure(figsize=(20,12))
        else:
            plt.figure()
        plt.plot(self.tvec,self.temp_cell, label="simulated cell temperature")
        plt.plot(self.tvec, self.temp_mud, label="simulated mud temperature")

        plt.plot(my_fluid_tvec,my_fluid_celltemps, '--', label="empirical cell temperature")
        plt.plot(my_fluid_tvec, my_fluid_temps, '--', label="empirical mud temperature")

        plt.legend()
        plt.ylabel("Temperature [K]")
        plt.xlabel("Time [s]")
        plt.title(f"{self.mud.name} with step input power of {self.input_power} W with Tf_final = {np.round(self.final_fluid_temperature,3)}")

    def __ss_model(self, t,x): # x = [Tc-Tamb; Tf-Tamb]
        
        Rj,Rfj,Rf = self.thermal_resistances
        Cc,Cf = self.thermal_capacitances


        A = np.array([[-(Rj+Rf)/(Cc*Rj*Rf), 1/(Rf*Cc)],
                [1/(Rf*Cf), - (1/(Rf*Cf) + 1/(Cf*Rfj))]])
        B = np.array([1/Cc, 0]).transpose() 
        C = np.array([0, 1])

        

        if self.controller_on:
            p = self.controller.feedback(x)
        
        else:
            if t > self.t_initial and t < self.t_stepfinal:
                p = self.input_power
            else:
                p = 0
        
        
        d = 0
        if p > 0:
            p = min(abs(p),self.parameters.max_heater)
        else:
            p = 0
        xdot = A @ x + B*p + d

        #print(x.shape, xdot)
        return xdot
