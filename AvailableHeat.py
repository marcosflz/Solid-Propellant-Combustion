# SCRIPT PARA OBTENER LA TEMPERATURA ADIABATICA DE LLAMA (METODO CALOR DISPONIBLE)

# En esta seccion se importa el modulo solidRocketData y la biblioteca NumPy.
import solidRocketData as srd  
import numpy as np

# Se define la funcion calorDisp la cual calcula la entalpia de formacion integrando el calor especifico mediante un metodo compuesto.
def calorDisp(h0,cp,t0,n):
    h = 10**3
    return n*(h0 + ((t0 - 298)/h)*(((cp(298) + cp(t0))/2) + sum([cp(298 + k*(t0 - 298)/h) for k in range(1,h-1)])))

# Se define la entalpia de formacion de los reactivos como la suma de sus componentes.
qReacts=sum([
calorDisp(srd.hfC12H22O11s,srd.cpC12H22O11,300,1),
calorDisp(srd.hfKNO3s,srd.cpKNO3,300,6.29)
])

# Se define el rango de temperaturas a analizar.
TList = np.arange(300,6000,1)

# Se definen las entalpias de formacion para el rango anterior.
qProd = [[
            calorDisp(srd.hfCO2g,srd.cpCO2,T,3.8),
            calorDisp(srd.hfCOg,srd.cpCO,T,5.21),
            calorDisp(srd.hfH2Og,srd.cpH2O,T,7.79),
            calorDisp(srd.hfH2g,srd.cpH2,T,3.07),
            calorDisp(srd.hfN2g,srd.cpN2,T,3.14),
            calorDisp(srd.hfK2CO3l,srd.cpK2CO3,T,3),
            calorDisp(srd.hfKOHl,srd.cpKOH,T,0.27)
        ]   for T in TList] 

# Se iteran las diferencias de entalpia en cada temperatura y se busca aquella diferencia que es minima definiendo su indice dentro del vector de datos
deltaQ = []
for i in range(0,len(TList)):
    deltaQ.append([TList[i],abs(qReacts-sum(qProd[i]))])

min_deltaQ = min(row[1] for row in deltaQ)
min_index = [row[1] for row in deltaQ].index(min_deltaQ)

# Se define la temperatura adiabatica de llama en base al criterio anterior.
TComb = deltaQ[min_index][0]

# Se definen otras propiedades de la reaccion y del gas formado por los productos
molecularMass_gas = ((3.8*44 + 5.21*28 + 7.79*18 + 3.07*2 + 3.14*28)/(3.8 + 5.21 + 7.79 + 3.07 + 3.14))*10**-3
cp_Molar_gas = (3.8*srd.cpCO2(TComb) + 5.21*srd.cpCO(TComb) + 7.79*srd.cpH2O(TComb) + 3.07*srd.cpH2(TComb) + 3.14*srd.cpN2(TComb))/(3.8 + 5.21 + 7.79 + 3.07 + 3.14)
cp_Mass_gas = cp_Molar_gas/molecularMass_gas
R_gas = 8.31446261815324/molecularMass_gas
gamma = cp_Mass_gas/(cp_Mass_gas - R_gas)


# Se sacan por pantalla los resultados.
print(
    ' T (K)              = ',deltaQ[min_index][0],'\n',
    'DeltaQ (J)         = ',deltaQ[min_index][1],'\n',
    'R (J/molK)         = ',R_gas,'\n',
    'gamma              = ',gamma,'\n',
    'Cp (J/molK)        = ',cp_Mass_gas,'\n',
    'M. Mass (kg/mol)   = ',molecularMass_gas,'\n',

)