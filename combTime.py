# SCRIPT PARA LA ESTIMACION DEL TIEMPO DE COMBUSTION

import math
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import scipy as sc
from scipy.optimize import fsolve
from scipy.interpolate import interp1d

#Datos de partida y referencia.

T1 = 1763
P3 = 91000
gamma = 1.2529712919314542 
R = 349.52459961214936 
h = 0.08
rho = 1667

# Areas reales de las gargantas de las toberas.
At_AS5 = 8.29044*10**-5
At_BN5 = 6.740926049286617*10**-5
At_MOC5 = 6.740926049286617*10**-5

regressionData = np.array([
    [68947.6,0.004],
    [344738,0.006],
    [689476,0.0075],
    [1.37895*10**6,0.0095],
    [2.41317*10**6,0.01125],
    [3.10264*10**6,0.013],
    [3.79212*10**6,0.012],
    [4.48159*10**6,0.013],
    [4.82633*10**6,0.0135],
    [5.17107*10**6,0.014],
    [5.86055*10**6,0.014],
    [5.86055*10**6,0.0155],
    [7.58424*10**6,0.015],
    [8.61845*10**6,0.0165],
    [9.30793*10**6,0.0175],
    [1.06869*10**7,0.017]
])

# Funcion para el calculo de ritmo de combustion del KNSU.
def regressionRate(P,a,n):
    return a*P**n
coefs, red = sc.optimize.curve_fit(regressionRate,regressionData[:,0],regressionData[:,1])

# Funcion que calcula de forma iterativa la combustion transitoria.
def combTime(At,R1i,R1f,s):

    R1 = np.arange(R1i,R1f+s,s)
    n = len(R1)

    Mi = abs(rho*np.pi*h*(R1f**2-R1**2))
    P1 = np.zeros(n)
    G = np.zeros(n)
    delta_m = np.zeros(n)
    delta_t = np.zeros(n)
    timeLine = np.zeros(n)

    C1 = 2*np.pi*rho*coefs[0]*h
    C2 = gamma*np.sqrt(((2/(gamma+1))**((gamma+1)/(gamma-1)))/(gamma*R*T1))

    for i in range(0,n-1):

        P1[i] = ((At/R1[i])*(C2/C1))**(1/(coefs[1]-1))
        G[i]  = C1*R1[i]*P1[i]**coefs[1]
        delta_m[i] = Mi[i] - Mi[i+1]
        delta_t[i] = delta_m[i]/G[i]
        

        timeLine = np.cumsum(delta_t)

    return  P1,G,delta_t,timeLine