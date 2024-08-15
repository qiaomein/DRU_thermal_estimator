%matplotlib widget

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

from parameters import *

Rc = heater_gap/(cell_area*k_ss)
Cf = 2000 * density_f * volume_f # estimate
Cc = cap_ss * m_ss
Rj = 3050 #0.03/(cell_area*0.3) # thermal resistance of the jacket
Rf = .3 #0.02/(cell_area*.4) # estimate equivalent thermal resistance relating to forced convection

P = 250 # 1000 Watts


# setup prelim transfer functions

tf2 = Rc*Cf*Cc
tf1 = Rc*Cf*(1/Rc + 1/Rj) + Cc*(1+Rc/Rf)
tf0 = (Rc/Rf + 1)*(1/Rc + 1/Rj) - 1/Rc
temp_fluid_TF = ct.TransferFunction([1],[tf2,tf1,tf0])
temp_cell_TF = ct.TransferFunction([Rc*Cf,Rc/Rf + 1],1) * temp_fluid_TF

# setup state space equations
A = np.matrix([[-1/Cc * (1/Rc + 1/Rj), 1/(Cf*Rc)],[1/(Cc*Rc), -1/Cf * (1/Rc + 1/Rf)]])
B = np.matrix([[1],[0]])
C = np.matrix([[1/Cc, 0],[0, 1/Cf]])

model = ct.StateSpace(A,B,C,0)

#%%
if __name__ == "__main__":
    print(f"values [Rc,Cc,Rf,Cf,Rj]: \n {[Rc,Cc,Rf, Cf, Rj]}")
    print("####################################################")
    print(f"A,B,C matrix: {A}, {B}, {C}\n")
    
    Tambient = 294
    N = 1000 # num of data points
    TN = 5000 #time span
    Tcell0 = 300
    Tfluid0 = 300
    Tref = 320 # mud temp should be here



    
    X0 = np.array([(Tcell0-Tambient)*Cc, (Tfluid0-Tambient)*Cf]) # qc, qf
    T = np.linspace(0,TN,N)
    
    U = P * np.concatenate((np.linspace(1,1,N-700), np.zeros(700)))
    yout, T, xout = lsim(model, U = U, T = T, X0=X0)
    yout += Tambient
    
    plt.figure()
    plt.plot(T,yout[:,0], label="cell temperature")
    plt.plot(T,yout[:,1], label="mud temperature")
    plt.legend()
    plt.ylabel("Temperature [K]")
    plt.xlabel("Time [s]")
    
    

    
    response = ct.step_response(P * model)
    plt.figure()
    response.plot()
    

