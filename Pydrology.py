#/usr/bin/python3
#Librería de métodos para modelación hidrológica SSIyAH-INA, 2022
import math
from sys import maxsize
from zlib import MAX_WBITS
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt
import os, glob

#0. Definición de funciones (métodos transversales)

def differentiate(list):
    dif=[0]*len(list)
    for i in range(1,len(list)):
        dif[i]=list[i]-list[i-1]
    return dif

#Computa Ecuación de Conservación
def waterBalance(Storage,Inflow,Outflow):
    Storage=Inflow-Outflow+Storage
    return Storage

#Computa EVR de acuerdo a las hipótesis de Thornthwaite
def computeEVR(P,EV0,Storage,MaxStorage):
    sigma=max(EV0-P,0)
    return(EV0+Storage*(1-math.exp(-sigma/MaxStorage))-sigma)

#Proratea Inflow
def apportion(Inflow,phi=0.1):
    return(Inflow*phi)

#Computa Runoff a partir del valor de Precipitación Efectiva Pe (acumulado en intervalo) y del almacenamiento a capacidad de campo o máximo almacenamiento
def curveNumberRunoff(NetRainfall,MaxStorage,Storage):
    return NetRainfall**2/(MaxStorage-Storage+NetRainfall)

#1. Objetos PQ: Componente de Producción de Escorrentía 

#1.A Reservorio de Detención 
class DetentionReservoir:
    """
    Reservorio de Detención. Vector pars de sólo parámetro: capacidad máxima de abstracción (MaxStorage). Vector de Condiciones Iniciales (InitialConditions): Storage. Condiciones de Borde (Boundaries): Inflow, EV, y Runoff). 
    """
    type='Detention Reservoir'
    def __init__(self,pars,InitialConditions=[0,0],Boundaries=[0,0],Proc='Abstraction'):
        self.MaxStorage=pars[0]
        self.Storage=InitialConditions[0]
        self.Runoff=InitialConditions[1]
        self.Inflow=Boundaries[0]
        self.EV=Boundaries[1]
        self.Proc=Proc
        self.dt=1
    def computeRunoff(self):
        if self.Proc == 'Abstraction':
            self.Runoff=max(0,self.Inflow-self.EV+self.Storage-self.MaxStorage)
        if self.Proc == 'CN_h0_continuous':
            self.Runoff=(max(self.Inflow*-self.EV,0))**2/(self.MaxStorage-self.Storage+self.Inflow-self.EV)
        self.Storage=waterBalance(self.Storage,self.Inflow,self.EV+self.Runoff)


#1.B Reservorio Lineal. 
class LinearReservoir:
    """
    Reservorio Lineal. Vector pars de un sólo parámetro: Tiempo de residencia (K). Vector de Condiciones Iniciales (InitialConditions): Storage y Outflow. Condiciones de Borde (Boundaries): Inflow y EV.
    """
    type='Linear Reservoir'
    def __init__(self,pars,InitialConditions=[0,0],Boundaries=[0,0],Proc='Agg',dt=1):
        self.K=pars[0]
        self.Storage=InitialConditions[0]
        self.Outflow=InitialConditions[1]
        self.Inflow=Boundaries[0] 
        self.EV=Boundaries[1]
        self.Proc=Proc
        if Proc == ('Agg' or 'API'):
            self.dt=1
        if Proc == 'Instant':
            self.dt=dt
    def computeOutFlow(self):
        if self.Proc == 'Agg':
            self.Outflow=(1/self.K)*(self.Storage+self.Inflow)
            self.Storage=waterBalance(self.Storage,self.Inflow,self.Outflow)
        if self.Proc == 'Instant':
            end=int(1/self.dt+1)    
            for t in range(1,end,1):
                self.Storage=waterBalance(self.Storage,self.Inflow*self.dt,self.Outflow*self.dt)
                self.Outflow=(1/self.K)*(self.Storage)
        if self.Proc == 'API':
            self.Storage=(1-1/self.K)*(self.Storage)+self.Inflow
            self.Outflow=(1/self.K)*self.Storage

class SCS_Reservoirs:
    type='Soil Conservation Service Model for Runoff Computation (Curve Number Method / Discrete Approach)'
    def __init__(self,pars,InitialConditions=[0,0,0],Boundaries=[0],Proc='Time Discrete Agg'):
        self.Abstraction=pars[0]
        self.MaxStorage=pars[1]
        self.Precipitation=Boundaries
        self.SoilStorage=InitialConditions[0]
        self.SurfaceStorage=[InitialConditions[1]]*len(self.Precipitation)
        self.Runoff=[InitialConditions[2]]*len(self.Precipitation)
        self.Infiltration=[0]*len(self.Precipitation)
        self.CumPrecip=[0]*len(self.Precipitation)
        self.NetRainfall=[0]*len(self.Precipitation) 
        self.Proc=Proc
        self.dt=1
    def computeAbstractionAndRunoff(self):
        for i in range(0,len(self.Precipitation)):
            if i == 0:
                  self.CumPrecip[i]=self.CumPrecip[i]+self.Precipitation[i]
            else:
                  self.CumPrecip[i]=self.CumPrecip[i-1]+self.Precipitation[i]
            if self.CumPrecip[i]-self.Abstraction > 0:
                   self.NetRainfall[i] = self.CumPrecip[i]-self.Abstraction
                   self.Runoff[i] = curveNumberRunoff(self.NetRainfall[i],self.MaxStorage,self.SoilStorage)
                   self.Infiltration[i]=self.NetRainfall[i]-self.Runoff[i]
            else:
                    self.NetRainfall[i] = 0
                    self.Runoff[i] = 0
            self.SurfaceStorage[i]=min(self.Abstraction,self.CumPrecip[i])
        self.Runoff=differentiate(self.Runoff)
        self.NetRainfall=differentiate(self.NetRainfall)
        self.Infiltration=differentiate(self.Infiltration)        

#2. Objetos PQ/QQ: Funciones de Distribución Temporal

#2.A Cascada de Reservorios Lineales (Discreta). Dos parámetros: Tiempo de Resdiencia (K) y Número de Reservorios (N)
class LinearReservoirCascade:
    """
    Cascada de Reservorios Lineales (Discreta). Vector pars de dos parámetros: Tiempo de Resdiencia (K) y Número de Reservorios (N). Vector de Condiciones Iniciales (InitialConditions): Si es un escalar genera una matriz de 2xN con valor constante igual al escalar. Condiciones de Borde (Boundaries): Inflow,Outflow (en este caso puede incluirse una matriz de caudales de 2XN)
    """
    type='Discrete Cascade of N Linear Reservoirs with Time K'
    def __init__(self,pars,Boundaries=[0,0],InitialConditions=[0],create='yes',Proc='Discretely Coincident',dt=1):
        self.K=pars[0]
        if not pars[1]:
            self.N=2
        else:
            self.N=pars[1]
        self.Inflow=Boundaries[0]   
        if  create == 'yes':
            self.Outflow=np.array([[InitialConditions[0]]*self.N]*2,dtype='float')
        else:
            self.Outflow=Boundaries[1]
        self.dt=dt
    def computeOutFlow(self):
        dt=self.dt
        k=self.K    
        c=math.exp(-dt/k)
        a=k/dt*(1-c)-c
        b=1-k/dt*(1-c)
        end=int(1/dt+1)
        for t in range(1,end,1):
            self.Outflow[1][0]=self.Inflow+(self.Outflow[0][0]-self.Inflow)*c
            if self.N > 1:
                for j in range(1,self.N,1):
                    self.Outflow[1][j]=c*self.Outflow[0][j]+a*self.Outflow[0][j-1]+b*self.Outflow[1][j-1]
            for j in range(0,self.N,1):
                self.Outflow[0][j]=self.Outflow[1][j]
    
if __name__ == "__main__":
    import sys

#3. Modelos PQ/QQ

# import Pydrology as Hydro 
# import matplotlib.pyplot as plt
# def TestRun(init=4,iter=100):
#     Cascade=Hydro.LinearReservoir_Cascade(pars,dt=0.01)
#     Cascade.ComputeOutflow()
#     vals=list()
#     vals.append(Cascade.Outflow[1][4])
#     Inflows=(init,0)
#     for t in range(1,iter):
#         if t > len(Inflows):
#             Cascade.Inflow=0
#         else:
#             Cascade.Inflow=Inflows[t-1]
#         Cascade.ComputeOutflow()
#         vals.append(Cascade.Outflow[1][4])
#     return(vals)


    
