# Solid-Propellant-Combustion
Python scripts for the design, analysis and prediction of the combustion of a solid propellant for rocket applications. 

The solidRocketData script is a "self-made" package in which all the data needed is included so that the Available Heat Method could be applied, it consists of different chemical and thermodynamic data for each reactant and product of a KNO3+Succrose combustion. The AvailableHeat script applies the self named method to this combustion to get the adiabatic combution temperature.

The combTime script uses external data of an experimental test to characterise the evolution of the burning rate with the pressure of the chamber. Taking into account all these data, a very simplified method can be written to simulate the evolution of the grain during the propellant combustion obtaining the transient propierties of the exhausted gases.


![alt text](https://github.com/marcosflz/Solid-Propellant-Combustion/blob/main/Images/Curva_Regresion_Interpolada.png)
![alt text](https://github.com/marcosflz/Solid-Propellant-Combustion/blob/main/Images/G_time.png)
![alt text](https://github.com/marcosflz/Solid-Propellant-Combustion/blob/main/Images/P1_time.png)
![alt text](https://github.com/marcosflz/Solid-Propellant-Combustion/blob/main/Images/T_time.png)
