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
        Tf = x[2]
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
        
        # find where the step input ends
        try:
            self.t_stepfinal = self.tvec[np.where(mud.heater_dutycycle > 0)[0][-1]] - self.mud.tvector[np.where(self.mud.heater_dutycycle > 0)[0][0]]
        except IndexError:
            print("NOT A STEP RESPONSE INPUT! Just lettin' ya know. :)")

        # parameters
        density_f = 1000 # water
        Rc = heater_gap/(cell_area*k_ss)

        Cf = 4000 * density_f * volume_f # estimate
        Cc = cap_ss * mass_ss 
        Ch = 1000
        Rc = 1
        Rj = 1.86 #0.03/(cell_area*0.3) # thermal resistance of the jacket
        Rf = 1.0428 #.1/(.6*cell_area) #0.02/(cell_area*.4) # estimate equivalent thermal resistance relating to forced convection
        Rfj = 1
        self.thermal_resistances = [Rc, Rj,Rfj,Rf]
        self.thermal_capacitances = [Ch, Cc, Cf]


        Rc,Rj,Rfj,Rf = self.thermal_resistances
        Ch,Cc,Cf = self.thermal_capacitances


        self.A = np.array([[-1/Ch * (1/Rc + 1/Rj), 1/(Rc*Ch), 0],
                [1/(Rc*Cc), -1/Cc * (1/Rc + 1/Rf), 1/(Rf*Cc)],
                [0, 1/(Rf*Cf), -1/Cf * (1/Rfj + 1/Rf)]])
        self.B = np.array([1/Ch, 0, 0]).transpose() 
        self.C = np.array([1, 0, 1]) # we can only sense Th and Tf

    def Tf_final(self, power):
        # final value theorem
        rc, rj,rfj,rf = self.thermal_resistances
        #cc,cf = self.thermal_capacitances

        bulk = rfj*rj/((rf+rj)*(rf+rfj))

        return power * bulk

    def stepResponse(self, x0 = None, powerPercentage = 1): # powerPercentage should be between 0-1
        

        assert powerPercentage > 0

        self.input_power = powerPercentage * self.parameters.max_heater


        # takes in true initial condition x0 = [tc,tf] 
        X0 = x0 - self.Tambient
        if x0 is None:
            X0 = np.zeros(3)
        sol = solve_ivp(self.__ss_model_step,(self.t_initial,self.t_final), X0, t_eval=self.tvec)

        yout = sol.y
        self.temp_mud = yout[2,:] + self.Tambient
        self.temp_cellwall = yout[1,:] + self.Tambient
        self.temp_cell = yout[0,:] + self.Tambient
        print(f"###############  {self.mud.name}   ################")
        print(f"Rc, Rj, Rfj, Rf: {self.thermal_resistances}")
        print(f"Ch, Cc, Cf: {self.thermal_capacitances}")
        print("###################################################")

        self.final_fluid_temperature = self.Tf_final(self.input_power) + self.Tambient

        return self.temp_cell, self.temp_mud
    
    def forcedResponse(self, x0 = None, pHistory = None):
  
        # takes in true initial condition x0 = [tc,tf] 
        X0 = x0 - self.Tambient
        if x0 is None:
            X0 = np.zeros(3)
        if pHistory is None:
            pHistory = self.mud.heater_dutycycle/100 * self.parameters.max_heater


        sol = solve_ivp(self.__ss_model,(self.t_initial,self.t_final), X0, t_eval=self.tvec)

        yout = sol.y
        self.temp_mud = yout[2,:] + self.Tambient
        self.temp_cellwall = yout[1,:] + self.Tambient
        self.temp_cell = yout[0,:] + self.Tambient
        print(f"###############  {self.mud.name}   ################")
        print(f"Rc, Rj, Rfj, Rf: {self.thermal_resistances}")
        print(f"Ch, Cc, Cf: {self.thermal_capacitances}")
        print("###################################################")

        self.final_fluid_temperature = self.temp_mud[-1]

        return self.temp_cell, self.temp_mud

    
    def plot(self, newfig = True, legend_on = True, stepResponse = False):
        t1 = self.mud.tvector[np.where(self.mud.heater_dutycycle > 0)[0][0]] # start when input is turned on
        if not stepResponse:
            t1 = 0
        t2 = self.mud.tvector[-1]
        self.plotting_indices = ( self.mud.tvector <= t2 ) & (self.mud.tvector >= t1)
        my_fluid_tvec = self.mud.tvector[self.plotting_indices] - t1

        my_fluid_temps = self.mud.temp_fluid[self.plotting_indices]
        my_fluid_celltemps = self.mud.temp_cell[self.plotting_indices]


        if newfig:
            plt.figure(figsize=(10,5))
        
        plt.subplot(2,1,1)
        plt.plot(self.tvec,self.temp_cell, label="simulated cell temperature")
        plt.plot(self.tvec,self.temp_cellwall, label="simulated cell wall temperature")
        plt.plot(self.tvec, self.temp_mud, label="simulated mud temperature")

        plt.plot(my_fluid_tvec,my_fluid_celltemps, '--', label="empirical cell temperature")
        plt.plot(my_fluid_tvec, my_fluid_temps, '--', label="empirical mud temperature")

        
        plt.ylabel("Temperature [K]")
        plt.xlabel("Time [s]")
        if stepResponse:
            plt.title(f"{self.mud.name} with step input power of {self.input_power} W with Tf_final = {np.round(self.final_fluid_temperature,3)} K")
        if legend_on:
            plt.legend(bbox_to_anchor=(1.05, 1.0))
        plt.tight_layout()

        plt.subplot(2,1,2)
        plt.plot(self.mud.tvector, self.mud.heater_dutycycle)

    def transfer_function(self):
        a = [list(l) for l in self.A]
        b = [[k] for k in self.B]
        c = list(self.C)

        return sig.ss2tf(a,b,c,0)
        
    def __ss_model(self, t,x): # x = [Tc-Tamb; Tf-Tamb]
        
        Rc,Rj,Rfj,Rf = self.thermal_resistances
        Ch,Cc,Cf = self.thermal_capacitances


        self.A = np.array([[-1/Ch * (1/Rc + 1/Rj), 1/(Rc*Ch), 0],
                [1/(Rc*Cc), -1/Cc * (1/Rc + 1/Rf), 1/(Rf*Cc)],
                [0, 1/(Rf*Cf), -1/Cf * (1/Rfj + 1/Rf)]])
        self.B = np.array([1/Ch, 0, 0]).transpose() 
        self.C = np.array([1, 0, 1]) # we can only sense Th and Tf

        power_history = self.mud.heater_dutycycle[t > self.mud.tvector]
        if len(power_history) == 0:
            p = 0
        else:

            p = self.mud.heater_dutycycle[t > self.mud.tvector][-1]/100 * self.parameters.max_heater
        
        
        d = 0
        if p > 0:
            p = min(abs(p),self.parameters.max_heater)
        else:
            p = 0
        xdot = self.A @ x + self.B*p + d

        #print(x.shape, xdot)
        return xdot


    def __ss_model_step(self, t,x): # x = [Tc-Tamb; Tf-Tamb]
        
        Rc,Rj,Rfj,Rf = self.thermal_resistances
        Ch,Cc,Cf = self.thermal_capacitances


        self.A = np.array([[-1/Ch * (1/Rc + 1/Rj), 1/(Rc*Ch), 0],
                [1/(Rc*Cc), -1/Cc * (1/Rc + 1/Rf), 1/(Rf*Cc)],
                [0, 1/(Rf*Cf), -1/Cf * (1/Rfj + 1/Rf)]])
        self.B = np.array([1/Ch, 0, 0]).transpose() 
        self.C = np.array([1, 0, 1]) # we can only sense Th and Tf

    
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
        xdot = self.A @ x + self.B*p + d

        #print(x.shape, xdot)
        return xdot
