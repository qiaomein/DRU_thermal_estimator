
# DRU Thermal Estimator

The goal of this project is to extract thermal properties of drilling muds while utilizing the existing hardware of the DRU Plus.

There are two main files: `dru_thermal_model.py` and `data_visualizer.py`.

## `data_visualizer.py`
This file handles the file parsing of data files from INSITE (the proprietary data logger) in csv format. There is the `Timestamp` class and the `MysteriousMud` class. 

The `Timestamp` class is strictly for taking a timestamp that is logged by INSITE and converting it into a time vector for plotting and data analysis. 

The `MysteriousMud` class reads in pertinent, raw data for relevant fields for a sample data. In this project, this includes dial readings, fluid temperature, cell temperature, heater duty cycle, and the time vector.

## `dru_thermal_model.py`

This file contains the mathematical model and all things related to the simulation of the thermal profile. 

The `ThermalModel` class contains the main 3rd order state space model along with methods for forced response fits. The system constant parameters can be set outside the class. The `ThermalController` class is an untested class that models the PID controller in the DRU Plus temperature conditioning phase. The original purpose of this was to emulate the controller and thermal dynamics together and perhaps design some kind of optimal controller. The `SystemParameters` class is a reference object that would contain theoretical values based on the geoemtry of the cell. Most of these values weren't used in the final model.
## Authors

- [@qiaomein](https://www.github.com/qiaomein)


## Installation

Simply install your latest version of Python 3 (this project was developed on Python 3.10.11) and install the following libraries

```bash
  pip install numpy
  pip install matplotlib
  pip install scipy
```
    
## Usage/Examples

```python
from dru_thermal_model import *
from data_visualizer import *

# read in test data
plt.figure(figsize=(10,5))
filename = "mud_pid1.csv"
sample_name = "BARAECD Fluid"
my_fluid = MysteriousMud(sample_name) # name the sample
my_fluid.importRawData(filename) # load in raw data

model = ThermalModel(my_fluid) # instantiate the model

# self.thermal_resistances = [Rc, Rj,Rfj,Rf]
# self.thermal_capacitances = [Ch, Cc, Cf]

# externally set the system parameters
model.thermal_resistances = [.13, 1.2, 2.6, .0018] # Rf
model.thermal_capacitances = [2000, 3000, 1000]
model.Tambient = 300 # set the ambient temperature

# set initial conditions
x0 = np.array([my_fluid.temp_cell[0],my_fluid.temp_fluid[0],my_fluid.temp_fluid[0]]) 

model.forcedResponse(x0) # run the simulation

# plot empirical data compared with simulation values
model.plot()
```


## Support

This is by no means a complete documentation of the libarry, but feel free to contact me for help/assistance at jackqiao2002@gmail.com

