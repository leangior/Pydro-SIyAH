#/usr/bin/python3
#Librería de métodos para modelación hidrológica SSIyAH-INA, 2022
import math
from typing import Optional, Union, List, Tuple
from sys import maxsize
from zlib import MAX_WBITS
from .pydrology_procedure_interface import PydrologyProcedureInterface
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt
import os, glob
import logging


#Funciones/Procedimientos (Procesos Hidrológicos, procesamiento de datos)
def testPlot(Inflow: Union[np.ndarray,List[float]],Outflow : Union[np.ndarray,List[float]]):
    """Genera una gráfica de prueba comparando 2 señales

     Args:
        Inflow : Union[np.ndarray,List[float]] 
            Hidrograma, lista de números (float)
        Outflow : Union[np.ndarray,List[float]]
            Hidrograma, lista de números (float)
    """
    plt.plot(Inflow,'r')
    plt.plot(Outflow,'b')
    plt.show()

def shiftLeft(array1d : Union[np.ndarray,List[float]],fill : float =0) -> np.ndarray:
    """Desplaza una serie hacia la izquierda (valor de índice en lista)

     Args:
        array1d : Union[np.ndarray,List[float]] 
            Serie numérica
        fill : float 
            valor de relleno para datos nulos
        
        Returns:
        Devuelve la lista original con el índice desplazado hacia la izquiera
    """
    i=0
    shift=np.array([0]*len(array1d),dtype='float')
    for value in range(1,len(array1d)):
        shift[i]=array1d[i+1]
        i=i+1
    shift[i]=fill
    return(shift)

def getDrivers(file : str,tCol: str ='t') -> pd.DataFrame:
    """Dummy para importar series temporales almacenadas en archivos CSV

     Args:
        file : str 
            ruta a archivo CSV
        tCol : str 
            nombre de columna con índices temporales (fecha/hora)
        
        Returns:
        Devuelve un dataframe python (serie temporal, indexada)
    """
    data=pd.read_csv(file)
    data[tCol]=pd.to_datetime(data[tCol],format='%Y-%m-%d')
    data.index=data[tCol]
    return(data)

def makeBoundaries(p : Union[List[float],np.ndarray] = [0],evp : Union[List[float],np.ndarray] =[0]) -> np.ndarray:
    """Dummy para generación de series de borde en modelos PQ operacionales (para cada polígono)

     Args: 
        p : [List[float],np.ndarray]  
            serie de datos de precipitación acumulada
        evp : [List[float],np.ndarray]  
            serie de datos de evapotranspiración potencial acumulada
        
        Returns:
        Devuelve un array 2d 
    """
    boundaries=np.array([[0]*2]*len(p),dtype='float')
    boundaries[:,0]=p
    boundaries[:,1]=evp
    return(boundaries)

def differentiate(lista_valores : List[float], asume_initial : bool =True) -> List[float]:
    """Diferencia la serie de valores
    
    Args:
        lista_valores : List[float] 
            lista de números (float)
        asume_initial : bool 
            indica si el primer elemento de la lista es el valor inicial, caso contrario asume que el valor inicial es 0
    
    Returns:
        Devuelve una lista con los valores diferenciados
    """
    if(asume_initial==True):
        dif=[lista_valores[0]]*len(lista_valores)
    else:
        dif=[0]*len(lista_valores)
    for i in range(1,len(lista_valores)):
        dif[i]=lista_valores[i]-lista_valores[i-1]
    return(dif)

def integrate(lista_valores : List[float], dt : float) -> float:
    """Integra por método del trapecio la serie de valores
    
    Args:
        lista_valores : List[float] 
            lista de números (float)
        dt : float 
            longitud de paso de cómputo
    
    Returns:
        Devuelve una lista con los valores integrados por método trapecio
    """
    integral=0
    for i in range(1,len(lista_valores)):
        integral=integral+(lista_valores[i]+lista_valores[i-1])*dt/2
    return(integral)

def triangularDistribution(pars : Union[List[float],float],distribution : str ='Symmetric',dt : float =0.01,shift : bool =True,approx: bool=True) -> np.ndarray:
    """Computa Hidrogramas Unitarios Triangulares (función respuesta a pulso unitario, función de transferencia) 
    
    Args:
        pars :  Union[List[float],float]
            lista de parámetros (tiempo al pico)
        distribution : str 
            tipo de distribución (simétrica, método SCS, asimétrica 'pbT')
    
    Returns:
        Devuelve un array1d con las ordenadas del HU
    """
    if isinstance(pars,(list)):
        T=pars[0]
    elif isinstance(pars,(float,int)):
        T=pars
    else:
        raise TypeError("pars must be a list, a float or an int")
    if distribution == 'Symmetric':
        tb=2*T
        peakValue=1/T
    elif distribution == 'SCS':
        tb=8/3*T
        peakValue=3/4*1/T
    elif distribution == 'pbT':
        tb=pars[0]
        peakValue=2/tb
    else:
        raise ValueError("distribution must be 'Symmetric', 'SCS' or 'pbT'")
    ndimu=int(round(tb/dt,0)+1)
    ndimU=int(round(tb,0)+1)
    u=np.array([0]*(ndimu),dtype='float')
    U=np.array([0]*(ndimU),dtype='float')
    i=0
    j=0    
    for t in np.array(list(range(0,ndimu))):
        if t*dt<T:
            if j==0:
                m=peakValue/T
            u[j]=m*t*dt
        else:
            if i==0:
                i=1
                m=peakValue/(T-tb)
            u[j]=peakValue+m*(t*dt-T)
        j=j+1 
    for j in range(1,ndimU):
        min=int((j-1)/dt)
        max=int(j/dt)+1
        U[j]=integrate(u[min:max],dt)
    if approx:
        U=U/sum(U)
    if shift:
        U=shiftLeft(U)         
        return(U[0:(len(U)-1)])
    else:
        return(U)

def gammaDistribution(n : float,k : float,dt : float =1,m : float =10,approx : bool = True,shift :bool =True) -> np.ndarray:
    """Computa Hidrogramas Unitarios (función respuesta a pulso unitario, función de transferencia) sobre la base de una función respuesta a impulso unitaria suponiendo n reservorios lineales en cascada con tiempo de residencia k (función de transferencia tipo gamma)
    
    Args:
        n : float
            cantidad de reservorios en la cascada (admite números reales)    
        k : float  
            tiempo de residencia
        m : 
            factor de tiempo de basse (longitud de base, expresada como 'mxn')
        approx: bool
            en caso de 'True', fuerza que la integral sea igual a 1
        shift: 
            en caso de 'True', desplaza a la izquierda a las ordenadas 

    Returns:
        Devuelve un array1d con las ordenadas del HU
    """
    T=int(m*n*k)
    u=np.array([0]*(int(T/dt)+1),dtype='float')
    U=np.array([0]*(int(T)+1),dtype='float')
    j=0
    for t in np.array(list(range(0,int(T/dt)+1))):
        u[j]=1/(k*math.gamma(n))*(t*dt/k)**(n-1)*math.exp(-t*dt/k)
        j=j+1
    for j in range(1,int(T)+1):
        min=int((j-1)/dt)
        max=int(j/dt)+1
        U[j]=integrate(u[min:max],dt)
    if approx:
        U=U/sum(U)
    if shift:
        U=shiftLeft(U)
    return(U)

def grXDistribution(T : float,distribution : str ='SH1') -> np.ndarray:    
    """Computa Hidrogramas Unitarios (función respuesta a pulso unitario, función de transferencia) propuestos por SEMAGREF (GR4J, GP)
    
    Args:
        T : float
            tiempo al pico    
        distribution : str  
            tipo de distribución 'SH1' corresponde a HU de 'flujo demorado', 'SH2' coresponde a HU de 'flujo rápido' 

    Returns:
        Devuelve un array1d con las ordenadas del HU
    """    
    if distribution == 'SH1':
        tb=T
        k=1
    elif distribution == 'SH2':
        tb=2*T
        k=1/2
    else:
        raise ValueError("Argumento distribution inválido")
    ndimU=int(tb)+2
    U=np.array([0]*(ndimU),dtype='float')
    for t in np.array(list(range(0,ndimU))):
        if t<T:
            U[t]=k*(t/T)**(5/2)
        else:
            if t>T and t<tb:
                 U[t]=1-k*(2-t/T)**(5/2)        
            else:
                if t>tb:
                    U[t]=1
    u=differentiate(U,True)
    return(shiftLeft(u))

def getPulseMatrix(inflows : Union[float,List[float],np.ndarray],u : np.ndarray) -> np.ndarray:
    """Computa matriz de convolución a partir de una lista o array1d de 'Inflows' (hidrograma de entrada) y sobre la base de una función de transferenciaa o HU 'u'    
    Args:
        Inflows : Union[float,List[float],np.ndarray]
            lista o array1d con hidrograma de entrada    
        u : np.ndarray  
            función de transferencia (HU)     

    Returns:
        Devuelve un array2d con la matriz de convolución 
    """    
    n=len(inflows)
    m=len(u)
    rows=n+m-1
    I=np.array([[0]*m]*rows,dtype='float')
    k=0
    for col in range(0,int(m)):
        for row in range(0,int(rows)):
            if col>row:
                I[row,col]=0
            else:
                if row>=n+k:
                    I[row,col]=0
                else:
                    I[row,col]=inflows[row-k]
        k=k+1
    return(I)

def waterBalance(Storage: float =0,Inflow : float =0,Outflow : float =0) -> float:
    """Computa la ecuación de conservación del volumen (balance hídrico)   
    Args:
        Storage: float
            Almacenamiento
        Inflow : float
            Suma de entradas (en paso de cálculo)   
        Outflow : float
            Suma de salidas  (en paso de cálculo)   

    Returns:
        Devuelve el valor inicial de almacenamiento (a fin de paso de cálculo) 
    """
    Storage=Inflow-Outflow+Storage
    return(Storage)

def computeEVR(P : float,EV0 : float ,Storage : float,MaxStorage : float) -> float:
    """Computa la evapotranspiración real de acuerdo a las hipótesis de Thornthwaite, siguiendo la ecuación formulada por Giordano (2019)   
    Args:
        Storage: float
            Almacenamiento a inicios de paso de cálculo
        MaxStorage: float
            Almacenamiento máximo (parámetro del modelo)
        EV0 : float
            Evapotranspiración potencial durante paso de cálculo 
        P : float
            Precipitación acumulada durante paso de cálculo   

    Returns:
        Devuelve el valor de la evapotranspiración real acumulada durante el paso de cálculo 
    """
    sigma=max(EV0-P,0)
    return(EV0+Storage*(1-math.exp(-sigma/MaxStorage))-sigma)

def apportion(Inflow : Union[float,List[float],np.ndarray],phi : float=0.1) -> Union[float,List[float],np.ndarray]:
    """Proratea un hidrograma   
    Args:
        Inflow: Union[float,List[float],np.ndarray]
            Hidorgrama de entrada
    Returns:
        Devuelve el hidrograma prorateado 
    """
    return(Inflow*phi)

def curveNumberRunoff(NetRainfall : float,MaxStorage : float,Storage : float) -> float:
    """Computa la escorrentía sobre la base de las hipótesis del método de SCS (Mockus, 1949)   
    
    Args:

        Storage: float
            Almacenamiento a inicios de paso de cálculo
        
        MaxStorage: float
            Almacenamiento máximo (parámetro del modelo)
        
        NetRainfall: 
            Precipitación neta durante el paso de cálculo
        
        Returns:
            Devuelve el valor de escorrentía producida durante el intervalo de cálculo
    """
    return NetRainfall**2/(MaxStorage-Storage+NetRainfall)


def SimonovKhristoforov(sim : np.array,obs : np.ndarray) -> np.ndarray: 
    """   
    Realiza correción de sesgo por simple updating (método propuesto por el Servicio Ruso)
    Args:

        sim: np.ndarray
            Serie simulada o sintética
        
        obs: np.ndarray
            Serie observada

        Returns:
            Devuelve serie simulada con correción de sesgo
    """
    uObs=np.mean(obs)
    uSim=np.mean(sim)
    df=np.array([[0]*2]*len(sim),dtype='float')
    df[:,0]=sim
    df[:,1]=obs
    df=pd.DataFrame(data=df)
    r=df.corr()[0][1]
    sSim=np.std(df[0])
    sObs=np.std(df[1])
    anomaly=sim-uSim
    return(uObs+r*sObs/sSim*anomaly)

#1. Proceso P-Q: Componentes de Función Producción de Escorrentía 

#1.A Reservorio de Retención 
class RetentionReservoir(PydrologyProcedureInterface):
    """
    Reservorio de Retención. Un sólo parámetro (capacidad máxima de almacenamiento). Condiciones de borde: lista con hidrograma/hietograma de entrada. Condición inicial       
    """
    type='Retention Reservoir'
    MaxStorage : float
    """Almacenamiento Máximo"""
    Inflow: np.ndarray
    """Serie de datos con hidrograma/hietograma de entrada (condición de borde)"""
    EVP: np.ndarray
    """Serie de datos de evapotranspiración potencial (condición de borde)"""
    Storage: np.ndarray
    """Almacenamiento a inicios de paso de cálculo (proceso computado)"""
    Runoff: np.ndarray
    """Escorrentía (proceso computado)"""
    EV: np.ndarray
    """Evapotranspiración real (proceso computado)"""
    Proc: str
    "Procedimiento: Abstraction o CN_h0_continuous"

    @property 
    def MaxStorage(self) -> float:
        return self.pars[0]
    
    def __init__(self,pars : List[float],InitialConditions : List[float] =[0],Boundaries : List[float] =[[0],[0]],Proc : str ='Abstraction'):
        """
            pars : list
                lista con el valor de almacenamiento máximo
            InitialConditions 
                lista con el valor de la condición inicial de almacenamiento
            Boundaries
                lista donde cada uno de los elementos es un vector columna con las condiciones de borde (Caudal Afluente/Precipitación y Evapotranspiración Potencial)
            Proc
                Procedimiento de cómputo para escorrentía. Admite 'Abstraction' o 'CN_h0_continuous'
            dt
                Longitud del paso de cómputo
        """
        super().__init__(pars,Boundaries,InitialConditions)
        self.Inflow=np.array(self.boundaries[0],dtype='float')
        self.EVP=np.array(self.boundaries[1],dtype='float')
        self.Storage=np.array([self.initial_conditions[0]]*(len(self.Inflow)+1),dtype='float')
        self.Runoff=np.array([0]*len(self.Inflow),dtype='float')
        self.EV=np.array([0]*len(self.Inflow),dtype='float')
        self.Proc=Proc
    def computeRunoff(self):
        for i in range(0,len(self.Inflow),1):
            self.EV[i]=computeEVR(self.Inflow[i],self.EVP[i],self.Storage[i],self.MaxStorage)
            if self.Proc == 'Abstraction':
                self.Runoff[i]=max(0,self.Inflow[i]-self.EV[i]+self.Storage[i]-self.MaxStorage)
            elif self.Proc == 'CN_h0_continuous':
                self.Runoff[i]=(max(self.Inflow[i]-self.EV[i],0))**2/(self.MaxStorage-self.Storage[i]+self.Inflow[i]-self.EV[i])
            else:
                raise ValueError("Argumento Proc inválido")
            self.Storage[i+1]=waterBalance(self.Storage[i],self.Inflow[i],self.EV[i]+self.Runoff[i])


#1.B Reservorio Lineal. ###REVISAR RUTINA TOMA LEN() OF UNSIZED OBJECT 
class LinearReservoir(PydrologyProcedureInterface):
    """
    Reservorio Lineal. Vector pars de un sólo parámetro: Tiempo de residencia (K). Vector de Condiciones Iniciales (InitialConditions): Storage, con el cual computa Outflow. Condiciones de Borde (Boundaries): Inflow y EV.
    """
    type='Linear Reservoir'
    K : float
    """ Constante de Recesión"""
    Inflow : np.ndarray 
    """ Caudal Alfuente"""
    EV: np.ndarray 
    """"Evpotranspiraciòn"""
    Storage : np.ndarray 
    """Almacenamiento"""
    Outflow: np.ndarray 
    """Caudal efluente"""

    @property
    def K(self) -> float:
        return self.pars[0]
    
    def __init__(self,pars : list,InitialConditions : list =[0],Boundaries : List[List[float]] =[[0],[0]],Proc : str ='Agg',dt : float=1):
        """
        pars : List[float]
            lista con el valor del coeficiente de recesión, expresado en dt unidades
        InitialConditions : list 
            lista con el valor de la condición inicial de almacenamiento
        Boundaries : List[List[float]]
            lista que contiene los vectores de las condiciones de borde (Caudal Afluente y Evapotranspiración Potencial)
        Proc : str
            Procedimiento de cómputo. Admite 'Agg' (Valor Agregado), 'API' o 'Instant' (Valor Instantáneo)
        dt : float
            Longitud del paso de cómputo
        """
        super().__init__(pars,Boundaries,InitialConditions)
        self.Inflow=np.array(self.boundaries[0],dtype='float') 
        self.EV=np.array(self.boundaries[1],dtype='float')
        self.Storage=np.array([self.initial_conditions[0]]*(len(self.Inflow)+1),dtype='float')
        self.Outflow=(1/self.K)*self.Storage
        self.Proc=Proc
        if Proc == ('Agg' or 'API'):
            self.dt=1
        elif Proc == 'Instant':
            self.dt=dt
        else:
            raise ValueError("Argumento Proc inválido")
    def computeOutFlow(self)-> None :
        for i in range (0,len(self.Inflow),1):
            if self.Proc == 'Agg':
                self.Outflow[i]=(1/self.K)*(self.Storage[i]+self.Inflow[i])
                self.Storage[i+1]=waterBalance(self.Storage[i],self.Inflow[i],self.Outflow[i])
            elif self.Proc == 'Instant':
                end=int(1/self.dt+1)
                Storage=self.Storage[i]
                Outflow=self.Outflow[i]    
                for t in range(1,end,1):
                    Storage=waterBalance(Storage,self.Inflow[i]*self.dt,Outflow*self.dt)
                    Outflow=(1/self.K)*(Storage)
                self.Storage[i+1]=Storage
                self.Outflow[i+1]=Outflow
            elif self.Proc == 'API':
                self.Storage[i+1]=(1-1/self.K)*self.Storage[i]+self.Inflow[i]
                self.Outflow[i+1]=(1/self.K)*self.Storage
            else:
                raise ValueError("Argumento Proc inválido")
                
class ProductionStoreGR4J(PydrologyProcedureInterface):
    """
    Reservorio de Producción de Escorrentía modelo GR4J
    """
    MaxSoilStorage : float
    """Almacenamiento Máximo en Perfil de Suelo"""
    Precipitation : np.ndarray
    """Precipitación (serie temporal)"""
    EVP : np.ndarray
    """Evapotranspiración potencial (serie temporal)"""
    SoilStorage : np.ndarray
    """Almacenamiento en Perfil de Suelo (serie temporal)"""
    NetEVP : np.ndarray
    """Evapotranspiración Potencial Neta (serie temporal)"""
    EVR : np.ndarray
    """Evapotranspiración Real (serie temporal)"""
    NetRainfall : np.ndarray
    """Precipitación Neta (serie temporal)"""
    Recharge : np.ndarray
    """Transferencia vertical: Pérdidas por recarga de flujo demorado (serie temporal)"""
    Infiltration : np.ndarray
    """Ingresos Netos al reservorio por Infiltración (serie temporal)"""
    Runoff : np.ndarray
    """Transferencia horizontal: Pérdidas por flujo directo (serie temporal)"""
    type='GR4J Runoff Production Store'

    @property
    def MaxSoilStorage(self) -> float:
        return self.pars[0]

    def __init__(self,pars : list,InitialConditions : list = [0],Boundaries : List[float] =[[0],[0]]):
        """
            pars : list
                lista con el valor de almacenamiento máximo
            InitialConditions : float
                lista con el valor de la condición inicial de almacenamiento
            Boundaries : List[float]
                lista de longitud 2 donde cada elemento es una lista que contiene los vectores de las condiciones de borde (Precipitación y Evapotranspiración Potencial)
        """
        super().__init__(pars,Boundaries,InitialConditions)
        self.Precipitation=np.array(self.boundaries[0],dtype='float')
        self.EVP=np.array(self.boundaries[1],dtype='float')
        self.SoilStorage=np.array([self.initial_conditions[0]]*(len(self.Precipitation)+1),dtype='float')
        self.NetEVP=np.array([0]*len(self.Precipitation),dtype='float')
        self.EVR=np.array([0]*len(self.Precipitation),dtype='float')
        self.NetRainfall=np.array([0]*len(self.Precipitation),dtype='float')
        self.Recharge=np.array([0]*len(self.Precipitation),dtype='float')
        self.Infiltration=np.array([0]*len(self.Precipitation),dtype='float')
        self.Runoff=np.array([0]*len(self.Precipitation),dtype='float')
    def computeOutFlow(self):
        for i in range(0,len(self.Precipitation)):
            self.NetRainfall[i]=max(0,self.Precipitation[i]-self.EVP[i])
            self.NetEVP[i]=max(0,self.EVP[i]-self.Precipitation[i])
            relativeMoisture=self.SoilStorage[i]/self.MaxSoilStorage
            ratio_netRainfall_maxStorage=self.NetRainfall[i]/self.MaxSoilStorage
            ratio_netEVP_maxStorage=self.NetEVP[i]/self.MaxSoilStorage
            self.Recharge[i]=self.MaxSoilStorage*(1-(relativeMoisture)**2)*np.tanh(ratio_netRainfall_maxStorage)/(1+relativeMoisture*np.tanh(ratio_netRainfall_maxStorage))
            self.EVR[i]=self.SoilStorage[i]*(2-relativeMoisture)*np.tanh(ratio_netEVP_maxStorage)/(1+(1-relativeMoisture)*np.tanh(ratio_netEVP_maxStorage))
            self.SoilStorage[i+1]=waterBalance(self.SoilStorage[i],self.Recharge[i],self.EVR[i])
            relativeMoisture=self.SoilStorage[i+1]/self.MaxSoilStorage
            self.Infiltration[i]=self.SoilStorage[i+1]*(1-(1+(4/9*relativeMoisture)**4)**(-1/4))
            self.SoilStorage[i+1]=waterBalance(self.SoilStorage[i+1],0,self.Infiltration[i])
            self.Runoff[i]=self.Infiltration[i]+self.NetRainfall[i]-self.Recharge[i]

class RoutingStoreGR4J(PydrologyProcedureInterface):
    """
    Reservorio de Propagación de Escorrentía modelo GR4J
    """
    type='GR4J Runoff Routing Store'
    MaxStorage : float
    """Almacenamiento Máximo"""
    waterExchange : float
    """Coeficiente de trasvase"""
    Inflow : np.ndarray
    """Recarga del reservorio (serie temporal)"""
    Leakages : np.ndarray
    """trasvase (serie temporal)"""
    Runoff: np.ndarray
    """Transferencia horizontal : pérdidas por flujo demorado"""
    Storage : np.ndarray
    """Almacenamiento"""

    @property
    def MaxStorage(self) -> float:
        return self.pars[0]
    
    @property
    def waterExchange(self) -> float:
        if len(self.pars)<2:
            return 0
        else:
            return self.pars[1]

    def __init__(self,pars : Union[List[float],List[Tuple[float,float]]],InitialConditions : list =[0],Boundaries : List[float] =[0]):
        """
            pars : Union[List[float],List[Tuple[float,float]]]
                lista con los valores de almacenamiento máximo (obligatorio) y del coeficiente de trasvase (opcional/puede omitirse, en cuyo caso se asigna el valor 0)
            InitialConditions 
                lista con el valor de la condición inicial de almacenamiento 
            Boundaries
                lista con el vector de la condición de borde (serie temporal de recarga de flujo demorado)
        """
        super().__init__(pars,Boundaries,InitialConditions)
        self.Inflow=np.array(self.boundaries[0],dtype='float')
        self.Leakages=np.array([0]*len(self.Inflow),dtype='float')
        self.Runoff=np.array([0]*len(self.Inflow),dtype='float')
        self.Storage=np.array([self.initial_conditions[0]]*(len(self.Inflow)+1),dtype='float')
    def computeOutFlow(self):
         for i in range(0,len(self.Inflow)):
            relativeMoisture=self.Storage[i]/self.MaxStorage
            self.Leakages[i]=self.waterExchange*relativeMoisture**(7/2)
            self.Storage[i+1]=max(0,self.Storage[i]+self.Inflow[i]+self.Leakages[i])
            relativeMoisture=self.Storage[i+1]/self.MaxStorage
            self.Runoff[i]=self.Storage[i+1]*(1-(1+relativeMoisture**4)**(-1/4))
            self.Storage[i+1]=waterBalance(self.Storage[i+1],0,self.Runoff[i])

class SCSReservoirs(PydrologyProcedureInterface):
    """
    Sistema de 2 reservorios de retención (intercepción/abstracción superficial y retención en perfil de suelo - i.e. capacidad de campo-), con función de cómputo de escorrentía siguiendo el método propuesto por el Soil Conservation Service. Lista pars de dos parámetros: Máximo Almacenamiento Superficial (Abstraction) y Máximo Almacenamiento por Retención en Perfil de Suelo (MaxStorage). Condiciones iniciales: Almacenamiento Superficial y Almacenamiento en Perfil de Suelo (lista de valores). Condiciones de Borde: Hietograma (lista de valores).
    """
    MaxSurfaceStorage : float
    """Almacenamiento Máximo en reservorio de retención"""
    MaxStorage : float
    """Almacenamiento Máximo en reservorio de producción"""
    Precipitation : np.ndarray
    """Precipitación (serie temporal)"""
    SurfaceStorage: np.ndarray
    """Almacenamiento en reservorio de retención (serie temporal)"""
    SoilStorage : np.ndarray
    """Almacenamiento en resservorio de producción (sserie temporal)"""
    Runoff : np.ndarray
    """Transferencia horizontal: escorrentía total (serie temporal)"""
    Infiltration : np.ndarray
    """Recarga de reservorio de producción (serie temporal)"""
    CumPrecip: np.ndarray
    "Precipitación acumulada durante el evento (serie temporal)"
    NetRainfall: np.ndarray
    """Precipitación neta (serie temporal)"""
    type='Soil Conservation Service Model for Runoff Computation (Curve Number Method / Discrete Approach)'
    
    @property
    def MaxSurfaceStorage(self) -> float:
        return self.pars[0]
    
    @property
    def MaxStorage(self) -> float:
        return self.pars[1]
    
    def __init__(self,pars : List[float],InitialConditions : List[float] =[0,0],Boundaries : List[float] =[0]):
        """
            pars : List[float]
                lista con los valores de almacenamiento máximo (reservorio de retención y reservorio de producción) y con el valor del coeficiente de recesión (flujo demorado) 
            InitialConditions : List[float]
                lista con los valores de la condición inicial de almacenamiento en cada reservorio
            Boundaries : List[float]
                lista con el vector de la condición de borde (serie temporal de precipitación)
        """
        super().__init__(pars,Boundaries,InitialConditions)
        self.Precipitation=np.array(self.boundaries,dtype='float')
        self.SurfaceStorage=np.array([self.initial_conditions[0]]*(len(self.Precipitation)+1),dtype='float')
        self.SoilStorage=np.array([self.initial_conditions[1]]*(len(self.Precipitation)+1),dtype='float')
        self.Runoff=np.array([0]*len(self.Precipitation),dtype='float')
        self.Infiltration=np.array([0]*len(self.Precipitation),dtype='float')
        self.CumPrecip=np.array([0]*len(self.Precipitation),dtype='float')
        self.NetRainfall=np.array([0]*len(self.Precipitation),dtype='float') 
    def computeAbstractionAndRunoff(self):
        Abstraction=self.MaxSurfaceStorage-self.SurfaceStorage[0]
        for i in range(0,len(self.Precipitation)):
            if i == 0:
                  self.CumPrecip[i]=self.CumPrecip[i]+self.Precipitation[i]
            else:
                  self.CumPrecip[i]=self.CumPrecip[i-1]+self.Precipitation[i]
            if self.CumPrecip[i]-Abstraction > 0:
                   self.NetRainfall[i] = self.CumPrecip[i]-Abstraction
                   self.Runoff[i] = curveNumberRunoff(self.NetRainfall[i],self.MaxStorage,self.SoilStorage[0])
                   self.Infiltration[i]=self.NetRainfall[i]-self.Runoff[i]
            else:
                    self.NetRainfall[i] = 0
                    self.Runoff[i] = 0
            self.SurfaceStorage[i+1]=min(self.SurfaceStorage[0]+Abstraction,self.SurfaceStorage[0]+self.CumPrecip[i])
        self.Runoff=differentiate(self.Runoff)
        self.NetRainfall=differentiate(self.NetRainfall)
        self.Infiltration=differentiate(self.Infiltration)
        for i in range(0,len(self.SoilStorage)-1):
            self.SoilStorage[i+1]=waterBalance(self.SoilStorage[i],self.Infiltration[i])        

class SCSReservoirsMod(PydrologyProcedureInterface):
    """
    Sistema de 2 reservorios de retención+detención (una capa de abstracción superficial/suelo y otra capa de retención/detención en resto perfil de suelo), con función de cómputo de escorrentía siguiendo el método propuesto por el Soil Conservation Service y añadiendo pérdida continua por flujo de base (primario). Lista pars de 3 parámetros: Máxima Abtracción por retención (Abstraction) y Máximo Almacenamiento por Retención+Detención en Perfil de Suelo (MaxStorage) y coefiente de pérdida K. Se añade pérdida continua. Condiciones iniciales: Almacenamiento Superficial y Almacenamiento en Perfil de Suelo (lista de valores). Condiciones de Borde: Hietograma (lista de valores).
    """
    MaxSurfaceStorage : float
    """Almacenamiento Máximo en reservorio de retención"""
    MaxStorage : float
    """Almacenamiento Máximo en reservorio de producción"""
    K : float
    """Coeficiente de recesión (autovalor dominante del sistema reservorio de produccción)"""
    Precipitation : np.ndarray
    """Precipitación (serie temporal)"""
    SurfaceStorage: np.ndarray
    """Almacenamiento en reservorio de retención (serie temporal)"""
    SoilStorage : np.ndarray
    """Almacenamiento en resservorio de producción (sserie temporal)"""
    Runoff : np.ndarray
    """Transferencia horizontal: flujo directo (serie temporal)"""
    Infiltration : np.ndarray
    """Recarga de reservorio de producción (serie temporal)"""
    CumPrecip: np.ndarray
    "Precipitación acumulada durante el evento (serie temporal)"
    NetRainfall: np.ndarray
    """Precipitación neta (serie temporal)"""
    BaseFlow: np.ndarray
    """Transferencia vertical: recarga de flujo demorado (serie temporal)"""
    type='Soil Conservation Service Model for Runoff Computation (Curve Number Method / Discrete Approach)'

    @property
    def MaxSurfaceStorage(self) -> float:
        return self.pars[0]
    
    @property
    def MaxStorage(self) -> float:
        return self.pars[1]
    
    @property
    def K(self) -> float:
        return self.pars[2]

    def __init__(self,pars : List[float],InitialConditions : List[float] =[0,0],Boundaries : List[float] =[0]):
        """
            pars : List[float]
                Lista con los valores de almacenamiento máximo (reservorio de retención y reservorio de producción) y con el valor del coeficiente de recesión (flujo demorado) 
            InitialConditions : Union[List[float],List[Tuple[float,float]]]
                lista de longitud 2 con el valor de la condición inicial de almacenamiento en cada reservorio
            Boundaries : List[float]
                lista con el vector de la condición de borde (serie temporal de precipitación)
        """
        super().__init__(pars,Boundaries,InitialConditions)
        self.Precipitation=np.array(self.boundaries,dtype='float')
        self.SurfaceStorage=np.array([self.initial_conditions[0]]*(len(self.Precipitation)+1),dtype='float')
        self.SoilStorage=np.array([self.initial_conditions[1]]*(len(self.Precipitation)+1),dtype='float')
        self.Runoff=np.array([0]*len(self.Precipitation),dtype='float')
        self.Infiltration=np.array([0]*len(self.Precipitation),dtype='float')
        self.CumPrecip=np.array([0]*len(self.Precipitation),dtype='float')
        self.NetRainfall=np.array([0]*len(self.Precipitation),dtype='float')
        self.BaseFlow=np.array([0]*len(self.Precipitation),dtype='float') 
    def computeAbstractionAndRunoff(self):
        Abstraction=self.MaxSurfaceStorage-self.SurfaceStorage[0]
        for i in range(0,len(self.Precipitation)):
            if i == 0:
                  self.CumPrecip[i]=self.CumPrecip[i]+self.Precipitation[i]
            else:
                  self.CumPrecip[i]=self.CumPrecip[i-1]+self.Precipitation[i]
            if self.CumPrecip[i]-Abstraction > 0:
                   self.NetRainfall[i] = self.CumPrecip[i]-Abstraction
                   self.Runoff[i] = curveNumberRunoff(self.NetRainfall[i],self.MaxStorage,self.SoilStorage[0])
                   self.Infiltration[i]=self.NetRainfall[i]-self.Runoff[i]
            else:
                    self.NetRainfall[i] = 0
                    self.Runoff[i] = 0
            self.SurfaceStorage[i+1]=min(self.SurfaceStorage[0]+Abstraction,self.SurfaceStorage[0]+self.CumPrecip[i])
        self.Runoff=differentiate(self.Runoff)
        self.NetRainfall=differentiate(self.NetRainfall)
        self.Infiltration=differentiate(self.Infiltration)
        for i in range(0,len(self.SoilStorage)-1):
            self.BaseFlow[i]=(1-self.K)*(self.SoilStorage[i]+self.Infiltration[i])
            self.SoilStorage[i+1]=waterBalance(self.SoilStorage[i],self.Infiltration[i],self.BaseFlow[i])        



#2. Proceso Q-Q: Componentes de Función Distribución de Escorrentía o Tránsito Hidrológico

#Cascada de Reservorios Lineales (Discreta). Dos parámetros: Tiempo de Resdiencia (K) y Número de Reservorios (N)
class LinearReservoirCascade(PydrologyProcedureInterface):
    """
    Cascada de Reservorios Lineales (Discreta). Lista pars de dos parámetros: Tiempo de Residencia (K) y Número de Reservorios (N). Vector de Condiciones Iniciales (InitialConditions): Si es un escalar (debe ingresarse como elemento de lista) genera una matriz de 2xN con valor constante igual al escalar, también puede ingresarse una matriz de 2XN que represente el caudal inicial en cada reservorio de la cascada. Condiciones de Borde (Boundaries): vector Inflow. 
    """
    K : float
    """Tiempo de residencia en reservorio"""
    N : float
    """Número de reservorios en cascada"""
    dt : float
    """Longitud del paso de cálculo"""
    type='Discrete Cascade of N Linear Reservoirs with Time K'

    @property
    def K(self) -> float:
        return self.pars[0]
    
    @property
    def N(self) -> float:
        if(len(self.pars)<2):
            return 2
        else:
            return self.pars[1]

    def __init__(self,pars : List[float],Boundaries : List[float] =[0],InitialConditions : List[float] =[0],dt=1):
        """
            pars : List[float]
                lista con los valores del tiempo de residencia K y de la cantidad discreta de N de reservorios lineales en cascada 
            InitialConditions : List[float]
                lista con el valor de caudal inicial y final para cada reservorio (puede brindarse un valor solamente, común a todos los reservorios, por defecto si se omite este es igual a 0)
            Boundaries : List[float]
                lista con el vector de la condición de borde (serie temporal de caudal afluente)
        """
        super().__init__(pars,Boundaries,InitialConditions)
        self.Inflow=np.array(self.boundaries)   
        if len(self.initial_conditions)==1:
            self.Cascade=np.array([[self.initial_conditions[0]]*self.N]*2,dtype='float')
            self.Outflow=np.array([self.initial_conditions[0]]*(len(self.boundaries)+1),dtype='float')
        else:
            if len(self.initial_conditions)!=2*self.N:
                raise NameError("Initial conditions list must have length 2*N")
            else:
                self.Cascade=np.array(self.initial_conditions,dtype='float').reshape(-1,2)
                self.Outflow=np.array([self.Cascade[0,1]]*(len(self.boundaries)+1),dtype='float')
        self.dt=dt
    def computeOutFlow(self):
        dt=self.dt
        k=self.K    
        c=math.exp(-dt/k)
        a=k/dt*(1-c)-c
        b=1-k/dt*(1-c)
        end=int(1/dt+1)
        for i in range(0,len(self.Inflow)):
            for n in range(1,end,1):
                self.Cascade[1][0]=self.Inflow[i]+(self.Cascade[0][0]-self.Inflow[i])*c
                if self.N > 1:
                    for j in range(1,self.N,1):
                        self.Cascade[1][j]=c*self.Cascade[0][j]+a*self.Cascade[0][j-1]+b*self.Cascade[1][j-1]
                for j in range(0,self.N,1):
                    self.Cascade[0][j]=self.Cascade[1][j]
            self.Outflow[i+1]=self.Cascade[0][j]

#Reservorio de enrutamiento (HIDROSAT)
class HIDROSATPowerLawReservoir(PydrologyProcedureInterface):
    """
    Reservorio potencial para tránsito de escorrentía total (modelo HIDROSAT). Se establece ley de potencia entre caudal (Q) y almacenamiento (W) en instante dado, de forma tal que Q=Q0.(W/W0)^(gamma). W0 representa el almacenamiento mínimo o de referencia para una situación en que la planicie aluvial se encuentre completamente anegada. Q0 es el caudal de referencia para dicha situación, mientras gamma es un parámetro de forma que tiene efecto sobre la capacidad de abstracción/retención de la red y sobre las puntas del hidrograma simulado (mayores o menores que un reservorio lineal para caso gamma!=1). El método de computación considera un único paso de cálculo, y estima un valor de descarga a fin de paso de cómputo, a partir de las condiciones de borde conocidas, pudiéndose subdividir en 1/dt subpasos de cómputo, de resolución dt. Se asume que las condiciones de borde están compuestas por una lista que contiene 4 valores: (a) la precipitación directa sobre el reservorio, (b) la evaporación desde la planicie a la atmósfera y (c.1) la escorrentía total que alcanza el reservorio al inicio del paso de cómputo y (c.2) la escorrentía total que alcanza el reservorio al final del paso de cómputo (pudiéndose informar sólo c.1, i.e. caso constante). La condición inicial está compuesta por una lista con 2 valores: (a) el almacenamiento inicial y (b) el área anegada inicial, expresada como fracción de área del sistema hídrico considerado.   
    """
    W0 : float
    """Almacenamiento de referencia en planicie aluvial"""
    Q0 : float
    """Caudal de referencia"""
    gamma: float
    """Factor de forma (abstracción) en relación Q(W)"""
    epsilon: float
    """Error máximo tolerable en Newthon-Raphson"""
    dt : float
    """Longitud del paso de cálculo"""
    type='Power Law Reservoir w/ vertical Losses'
    
    @property
    def W0(self) -> float:
        return self.pars[0]
    
    @property
    def Q0(self) -> float:
        return self.pars[1]
    
    @property
    def gamma(self) -> float:
        return self.pars[2]
    
    @property
    def maxFlooded(self) -> float:
        if(len(self.pars)<4):
            return 1
        else:
            return self.pars[3]

    @property
    def epsilon(self) -> float:
        if(len(self.pars)<5):
            return 0.0001
        else:
            return self.pars[4]
    
    def __init__(self,pars : List[float],Boundaries : List[float] =[0,0,0,0],InitialConditions : List[float] =[0,0],dt=1):
        """
            pars : List[float]
                lista con los valores del almacenamiento de referencia W0, el caudal de referencia Q0, el factor de forma y opcionalmente el valor máximo de área inundada y el valor de epsilon 
            InitialConditions : List[float]
                lista con los valores de la condición inicial de almacenamiento el reservorio y el área anegada en la planicie aluvial
            Boundaries : List[float]
                lista con las condiciones de borde (precipitación directa, Evaporación y la escorrentía afluente a inicios y final de paso de cómputo
        """
        super().__init__(pars,Boundaries,InitialConditions)
        self.Precipitation=np.array(self.boundaries[0],dtype='float')
        self.EV=np.array(self.boundaries[1],dtype='float')
        if len(self.boundaries) < 4:
            self.Inflow=self.boundaries[2]
            self.Inflow=np.array([self.Inflow]*2,dtype='float')
        elif len(self.boundaries) == 4:
            self.Inflow=[self.boundaries[2],self.boundaries[3]]
            self.Inflow=np.array(self.Inflow,dtype='float')
        else:
            raise NameError("Boundaries must be a List object of length 3 or 4 (direct precipitation, evaporation, initial/final inflows of time step)")
        self.Storage=np.array([self.initial_conditions[0]]*2,dtype='float')
        self.Flooded=np.array([self.initial_conditions[1]]*2,dtype='float')
        self.Outflow=self.Q0*(self.initial_conditions[0]/self.W0)**(self.gamma)
        self.Outflow=np.array([self.Outflow]*len(self.Inflow),dtype='float')
        self.dt=dt

    def estimateDischarge(self):
        self.U=(self.Inflow[0]+self.Inflow[1]-self.Outflow[0])*self.dt/2+(self.Precipitation-self.EV)*self.dt+self.Storage[0]
        self.U=max(0,self.U)
        qPred=self.Q0*(self.U/self.W0)**(self.gamma)
        if(qPred>0):
            self.f=self.W0*(qPred/self.Q0)**(1/self.gamma)+qPred/2*self.dt-self.U
            self.df=self.W0/(self.Q0*self.gamma)*(qPred/self.Q0)**(1/self.gamma-1)+self.dt/2
        else:
            self.f=0
            self.df=self.dt/2
        self.Outflow[1]=max(0,qPred-self.f/self.df)
        while abs(qPred-self.Outflow[1])>self.epsilon:
            qPred=self.Outflow[1]
            if(qPred>0):
                self.f=self.W0*(qPred/self.Q0)**(1/self.gamma)+qPred/2*self.dt-self.U
                self.df=self.W0/(self.Q0*self.gamma)*(qPred/self.Q0)**(1/self.gamma-1)+self.dt/2
            else:
                self.f=0
                self.df=self.dt/2
            self.Outflow[1]=max(0,qPred-self.f/self.df)
        
    def computeOutFlow(self):
        if len(self.Inflow)!=2:
            raise NameError("Inflow array must have lenght 2")
        for t in range(1,round(1/self.dt+1)):
            self.estimateDischarge()
            self.Outflow[0]=self.Outflow[1]
            self.Storage[0]=self.W0*(float(self.Outflow[0])/self.Q0)**(1/self.gamma)
        self.Storage[0]=self.initial_conditions[0]
        self.Outflow[0]=self.Q0*(self.initial_conditions[0]/self.W0)**(self.gamma)
        self.Storage[1]=self.W0*(float(self.Outflow[1])/self.Q0)**(1/self.gamma)
        self.Flooded[1]=self.maxFlooded*min(1,(float(self.Storage[1])/self.W0)**(self.gamma))


#Canal Muskingum 
# EN DESARROLLO (MUSKINGUM y CUNGE) --> VER RESTRICCIONES NUMÉRICAS y SI CONSIDERAR CURVAS HQ y BH COMO PARAMETROS DEL METODO. POR AHORA FINALIZADO MUSKINGUM CLÁSICO. CUNGE DEBE APOYARSE SOBRE EL MISMO, MODIFICANDO PARS K y X
class MuskingumChannel(PydrologyProcedureInterface):
    """
    Método de tránsito hidrológico de la Oficina del río Muskingum. Lista pars de dos parámetros: Tiempo de Tránsito (K) y Factor de forma (X). Condiciones Iniciales (InitialConditions): lista con array de condiciones iniciales o valor escalar constante. Condiciones de borde: lista con hidrograma en nodo superior de tramo. A fin de mantener condiciones de estabilidad numérica en la propagación (conservar volumen), sobre la base de la restricción 2KX<=dt<=2K(1-X) (Chin,2000) y como dt viene fijo por la condición de borde (e.g. por defecto 'una unidad') y además se pretende respetar el valor de K, se propone incrementar la resolución espacial dividiendo el tramo en N' subtramos de igual longitud, con tiempo de residencia mínimo T'=K/N', para el caso dt<2KX (frecuencia de muestreo demasiado alta). Así para obtener el valor N', se aplica el criterio de Chin estableciendo que el valor crítico debe satisfacer dt=uT', específicamente con u=2X y T' = K/N'--> N'=2KX/dt. Al mismo tiempo si dt>2K(1-X) (frecuencia de muestreo demasiado baja), el paso de cálculo se subdivide en M subpasos de longitud dT=2K(1-X) de forma tal que dT/dt=dv y M=dt/dv. El atributo self.tau resultante especifica el subpaso de cálculo (siendo self.M la cantidad de subintervalos utilizados) y self.N la cantidad de subtramos. 
    """ 
    K : float
    """Tiempo de tránsito, parámetro del modelo"""
    X : float
    """Factor de forma, parámetro del modelo"""
    dt : float
    """Longitud del paso de cálculo"""
    Inflow: np.ndarray
    """"Hidrogama de ingresos al tramo (serie temporal)"""
    Outflow: np.ndarray 
    """Hidrograma de descargas del tramo (serie temporal)"""
    N : float
    """Cantidad de subtramos en tramo"""
    M : float
    """Cantidad de subpasos de cálculo en paso"""
    tau : float
    """Longitud de subpaso de cómputo"""
    type='Muskingum Channel'

    @property
    def K(self) -> float:
        return self.pars[0]
    
    @property
    def N(self) -> float:
        return self.pars[1]

    def __init__(
        self,
        pars : List[float],
        Boundaries : List[float] = [0],
        InitialConditions : List[float] =[0],
        dt : float = 1):
        """
            pars : List[float]
                lista con los valores del tiempo de tránsito (K) y del factor de forma (X) 
            InitialConditions : List [float]
                lista con los valores de la condición inicial de almacenamiento en tramo 
            Boundaries : List[float]
                lista con el hidrograma de entrada, de resolución dt
            dt : float
                resolución del hidrograma de entrada, por defecto se establece en la unidad
        """
        super().__init__(pars,Boundaries,InitialConditions)
        self.dt=dt
        self.lowerbound=2*self.K*self.X
        self.upperbound=2*self.K*(1-self.X)
        self.Inflow=np.array(self.boundaries,dtype='float')
        self.Outflow=np.array([0]*len(self.Inflow),dtype='float')
        self.InitialConditions=np.array(self.initial_conditions,dtype='float')
        self.N=1
        self.tau=self.dt
        if self.dt > self.upperbound:
            self.tau=self.upperbound
        else:
            if self.dt < self.lowerbound:
               self.N=round(self.lowerbound/self.dt) 
        self.M=round(self.dt/self.tau)
        if self.X > 1/2:
            raise NameError('X must be between 0 and 1/2')
        if len(InitialConditions) == 1:
            if InitialConditions[0] == 0:
                self.InitialConditions=np.array([[0]*(self.N+1)]*2,dtype='float')
            else:    
                self.InitialConditions=np.array([[InitialConditions[0]]*(self.N+1)]*2,dtype='float')
        if len(self.InitialConditions[0]) < self.N:
            raise NameError('Matrix of Initial Conditions must have'+str(self.N+1)+'cols as it have'+str(self.N)+'subreaches')        
        self.Outflow[0]=self.InitialConditions[1][self.N]
    def computeOutFlow(self):
        K=self.K/self.N
        X=self.X
        tau=self.tau
        D=(2*K*(1-X)+tau)    
        C0=(tau+2*K*X)/D
        C1=(tau-2*K*X)/D
        C2=(2*K*(1-X)-tau)/D
        for i in range(0,len(self.Inflow)-1,1):
            self.InitialConditions[0][0]=self.Inflow[i]
            self.InitialConditions[1][0]=self.Inflow[i+1]
            for j in range(1,self.N+1,1):
                for t in range(0,self.M,1):
                    self.InitialConditions[1][j]=C0*self.InitialConditions[0][j-1]+C1*self.InitialConditions[1][j-1]+C2*self.InitialConditions[0][j]
                    self.InitialConditions[0][j]=self.InitialConditions[1][j]
            self.Outflow[i+1]=max(self.InitialConditions[1][self.N],0)    
            

#Tránsito Lineal con funciones de transferencia. Por defecto, se asume una distrinución gamma con parámetros n (número de reservorios) y k (tiempo de residencia). Asimismo, se considera n=2, de modo tal que tp=k (el tiempo al pico es igual al tiempo de residencia) 
class LinearChannel(PydrologyProcedureInterface):
    """
    Método de tránsito hidrológico implementado sobre la base de teoría de sistemas lineales. Así, considera al tránsito de energía, materia o información como un proceso lineal desde un nodo superior hacia un nodo inferior. Específicamente, sea I=[I1,I2,...,IN] el vector de pulsos generados por el borde superior y U=[U1,U2,..,UM] una función de distribución que representa el prorateo de un pulso unitario durante el tránsito desde un nodo superior (borde) hacia un nodo inferior (salida), el sistema opera aplicando las propiedades de proporcionalidad y aditividad, de manera tal que es posible propagar cada pulso a partir de U y luego mediante la suma de estos prorateos obtener el aporte de este tránsito sobre el nodo inferior (convolución).
    """
    Inflow : np.ndarray
    """Hidrograma de entrada (serie temporal)"""
    dt : float
    """Longitud de paso de cálculo"""
    Proc: str
    """Tipo de Procedimiento. Admite 'Nash' (cascada de reservorios lineales, debe proveerse lista de pars k y n) y 'UH' (Hidrograma Unitario, debe proveerse lista con array conteniendo ordenadas de UH a paso dt)"""
    type='Single Linear Channel'
    def __init__(self,pars : List[float],Boundaries : List[float] =[0],Proc : str='Nash',dt : float =1):
       """
            pars : List[List[float]]
                Lista de flotantes con los valores del tiempo de residencia (k) y número de reservorios (n), en caso que Proc='Nash', o lista, tupla o array con ordenadas de Hidrograma Unitario, en caso que  Proc='UH'
            InitialConditions : List [float]
                Lista  con el valor de la condición inicial de almacenamiento en tramo 
            Boundaries : List[float]
                Lista con el hidrograma de entrada, de resolución dt
            dt : float
                Resolución del hidrograma de entrada, por defecto se establece en la unidad
            Proc: str
                Procedimiento para transferencia: 'Nash' (cascada de Nash, debe proveerse k y n) y 'UH' (Hidrograma Unitario, lista, tupla o array con valoress de ordenadass)

       """
       super().__init__(pars,Boundaries)
       self.routingProc=Proc
       self.pars=np.array(self.pars,dtype='float')
       self.Inflow=np.array(self.boundaries,dtype='float')
       self.dt=dt
       if self.routingProc == 'Nash':
            self.k=self.pars[0]
            self.n=self.pars[1]
            self.u=gammaDistribution(self.n,self.k,self.dt)
       elif self.routingProc == 'UH':
            self.u=self.pars
       else:
            raise ValueError("Argumento Proc inválido. Debe ser 'Nash' o 'UH'")
       self.Outflow=np.array([[0]]*(len(self.Inflow)+len(self.u)-1))
    def computeOutFlow(self):
        I=getPulseMatrix(self.Inflow,self.u)
        self.Outflow=np.dot(I,self.u)

class LinearNet(PydrologyProcedureInterface):
    """
    Método de tránsito hidrológico implementado sobre la base de teoría de sistemas lineales. Así, considera al tránsito de energía, materia o información como un proceso lineal desde N nodos superiores hacia un nodo inferior. Específicamente, sea I=[I1,I2,...,IN] un vector de pulsos generados por un borde y U=[U1,U2,..,UM] una función de distribución que representa el prorateo de un pulso unitario durante el tránsito desde un nodo superior (borde) hacia un nodo inferior (salida), aplicando las propiedades de proporcionalidad y aditividad es posible propagar mediante convolución cada pulso de cada hidrograma a partir de su respectiva función de distribución U y luego mediante la suma de las propagaciones obtenerse el aporte de este tránsito sobre el nodo inferior. Numéricamente el sistema se representa como una transformación matricial (matriz de pulsos*u=vector de aportes). Consecuentemente, el tránsito se realiza para cada borde y la suma total de estos tránsitos constituye la señal transitada sobre el nodo inferior.  Condiciones de borde: Lista con los hidrogramas en nodos superiores del tramo. Parámetros: función de distribución (proc='UH') o tiempo de residencia (k) y número de reservorios (n), si se desea utilizar el método de hidrograma unitario de Nash (proc='Nash'). Pars es una lista en donde la información necesaria para cada nodo se presenta por fila (parámetros de nodo). El parámetro dt refiere a la longitud de paso de cálculo para el método de integración, siendo dt=1 la resolución nativa de los hidrogramas de entrada provistos. Importante, las funciones de transferencia deben tener la misma cantidad de ordenadas (igual dimensión). 
    """
    pars : np.ndarray
    """Matriz con parámetros de propagación (j-vectores fila)"""
    Inflows : np.ndarray
    """Matriz con hidrogramas de entrada (j-vectores columna)"""
    Proc : str
    """Procedimiento, admite 'Nash' y 'UH'"""
    dt : float
    """Longitud del paso de cómputo"""
    type='Linear Routing System. System of Linear Channels'
    def __init__(self,pars : List[List[float]],Boundaries : List[List[float]] ,Proc : str = 'Nash',dt : float =1):
        """
            pars : List[List[float]]
                Lista  con los valores del tiempo de residencia (k) y número de reservorios (n), en caso que Proc='Nash', o con ordenadas de cada Hidrograma Unitario, en caso que  Proc='UH', para cada nodo de entrada
            Boundaries : List[List[float]]
                Lista con los hidrogramas de cada j-ésimo nodo de entrada 
            dt : float
                Resolución del hidrograma de entrada, por defecto se establece en la unidad
            Proc: str
                Procedimiento para transferencia: 'Nash' (cascada de Nash, debe proveerse k y n) o 'UH' (Hidrogramas Unitarios, array con j-vectores fila con valores de ordenadas)

        """
        super().__init__(pars,Boundaries)
        self.routingProc=Proc
        if not (self.routingProc=='Nash' or self.routingProc=='UH'):
           raise ValueError("Argumento Proc inválido. Debe ser 'Nash' o 'UH'")
        self.pars=np.array(self.pars,dtype='float')
        self.Inflows=np.array(self.boundaries,dtype='float')
        self.dt=dt
    def computeOutflow(self):
        for j in range(0,len(self.Inflows)):
            linear = LinearChannel(
                pars = self.pars[j],
                Boundaries = self.Inflows[j],
                dt = self.dt,
                Proc = self.Proc)
            linear.computeOutFlow()
            if j == 0:
                self.Outflow = linear.Outflow
            else:
                nrows = max(len(self.Outflow), len(linear.Outflow))
                f = np.zeros((nrows))
                if len(self.Outflow) > len(linear.Outflow):
                   f[0:len(linear.Outflow)] = linear.Outflow[0:len(linear.Outflow)]
                   self.Outflow=self.Outflow + f
                else:
                   f[0:len(self.Outflow)] = self.Outflow[0:len(self.Outflow)]
                   self.Outflow = f + linear.Outflow 

#Método de Clark
class ClarkSystem(PydrologyProcedureInterface):
    """
    Método de tránsito hidrológico implementado sobre la base de la propuesta de Clark. Requiere un UH () obtenido sobre la base de una curva TAC (análisis de MDET) y un tiempo de residencia (tr). El método es una extensión del enfoque modelístico, asumiendo la cuenca como sistema lineal y sobre la base de la convolución de los pulsos de escorrentía, utilizando una función de transferencia'tiempo-area' (TAC). Así, primeramente se obtiene un hidrograma por convolución (escorrentía/TAC) y este se transita sobre un reservorio lineal a fin de incluir el efecto de almacenamiento en la red de drenaje.  
    """
    Inflow : np.ndarray
    """Hidrograma de entrada (serie temporal)"""
    dt : float
    """Longitud de paso de cálculo"""
    type='Clark Routing System'
    
    @property
    def u(self) -> float:
        if self.pars[0] is not None:
            return self.pars[0]
        else:
            raise Exception("Ordinates of Time Area curve must be provided. Check pars[0]")

    @property
    def k(self) -> float:
        if self.pars[1] is not None:
            return self.pars[1]
        else:
            raise Exception("Residence time must ne provided. Check pars[1]")
        
    def __init__(self,pars : List[float],Boundaries : List[float] =[[0],[0]],InitialConditions : float = [0],Proc : str = 'Clark',dt=1):
       """
        pars : List[List[float]]
            Lista que contiene las ordenadas de la curva 'tiempo area' (lista de flotantes) y el valor del tiempo de residencia en el reservorio lineal (k)
        InitialConditions : List [float]
            Lista  con el valor de la condición inicial de almacenamiento en el reservorio lineal
        Boundaries : List[float]
            Lista con el hidrograma de entrada (resolución dt) y opcionalmente una lista con valores de pérdida (leakages)  
        dt : float
             Resolución del hidrograma de entrada, por defecto se establece en la unidad
       """
       super().__init__(pars,Boundaries,InitialConditions)
       self.routingProc=Proc
       self.Inflow=np.array(self.boundaries[0],dtype='float')
       if len(self.boundaries)>1:
            if self.boundaries[1]==len(self.boundaries[0]):
                    self.Leakages=np.array(self.boundaries[1],dtype='float')
            elif len(self.boundaries[1])<len(self.boundaries[0]):
                    dif=len(self.boundaries[0])-len(self.boundaries[1])
                    self.Leakages=np.append(self.boundaries[1],np.zeros(dif))
            else:
                dif=len(self.boundaries[1])-len(self.boundaries[0])
                self.Inflow=np.append(self.boundaries[0],np.zeros(dif))
                self.Leakages=self.boundaries[1]
       else:
            self.Leakages=np.zeros(len(self.boundaries[0]))
       self.dt=dt
       
    def executeRun(self):
       self.TimeAreaConvolution=LinearChannel(pars=self.u,Boundaries=self.Inflow,Proc='UH',dt=self.dt)
       self.TimeAreaConvolution.computeOutFlow()
       self.RoutingReservoir=LinearReservoir(pars=[self.k],InitialConditions=self.initial_conditions,Boundaries=[self.TimeAreaConvolution.Outflow,self.Leakages],Proc='Instant',dt=self.dt)
       self.RoutingReservoir.computeOutFlow()
       self.Q=self.RoutingReservoir.Outflow

class LagAndRoute(PydrologyProcedureInterface):
    """
    Método de tránsito hidrológico implementado sobre la base de la propuesta de Clark. Requiere un UH () obtenido sobre la base de una curva TAC (análisis de MDET) y un tiempo de residencia (tr). El método es una extensión del enfoque modelístico, asumiendo la cuenca como sistema lineal y sobre la base de la convolución de los pulsos de escorrentía, utilizando una función de transferencia'tiempo-area' (TAC). Así, primeramente se obtiene un hidrograma por convolución (escorrentía/TAC) y este se transita sobre un reservorio lineal a fin de incluir el efecto de almacenamiento en la red de drenaje.  
    """
    Inflow : np.ndarray
    """Hidrograma de entrada (serie temporal)"""
    dt : float
    """Longitud de paso de cálculo"""
    type='Lag and route method'
    
    @property
    def lag(self) -> float:
        if self.pars[0] is not None:
            return self.pars[0]
        else:
            raise Exception("lag time must be provided. Check pars[0]")

    @property
    def k(self) -> float:
        if self.pars[1] is not None:
            return self.pars[1]
        else:
            raise Exception("Residence time (k) must ne provided, for computation you may take attenuation index = -1/ln(k). Check pars[1]")
    
    @property
    def n(self) -> float:
        if len(self.pars)>2:
            print("assuming Cascade of Linear Reservoirs Routing")
            return self.pars[2]
        else:
            print("assuming Linear Reservoir Routing")
            return 1
        
    def __init__(self,pars : List[float],Boundaries : List[float] =[[0],[0]],InitialConditions : float = [0],Proc : str = 'Lag and Route',dt=1):
       """
        pars : List[List[float]]
            Lista que contiene el tiempo de retardo (celeridad), el tiempo de residencia (memoria) y opcionalmente el número de reservorios (lineales) utilizados en la propagación
        InitialConditions : List [float]
            Lista  con el valor de la condición inicial de almacenamiento en el sistema de reservorios (por defecto se cosnsideran vacíos)
        Boundaries : List[float]
            Lista con el hidrograma de entrada (resolución dt) y opcionalmente una lista con valores de pérdida (leakages)  
        dt : float
             Resolución del hidrograma de entrada, por defecto se establece en la unidad
       """
       super().__init__(pars,Boundaries,InitialConditions)
       self.routingProc=Proc
       self.Inflow=np.array(self.boundaries[0],dtype='float')
       if len(self.boundaries)>1:
            if self.boundaries[1]==len(self.boundaries[0]):
                    self.Leakages=np.array(self.boundaries[1],dtype='float')
            elif len(self.boundaries[1])<len(self.boundaries[0]):
                    dif=len(self.boundaries[0])-len(self.boundaries[1])
                    self.Leakages=np.append(self.boundaries[1],np.zeros(dif))
            else:
                dif=len(self.boundaries[1])-len(self.boundaries[0])
                self.Inflow=np.append(self.boundaries[0],np.zeros(dif))
                self.Leakages=self.boundaries[1]
       else:
            self.Leakages=np.zeros(len(self.boundaries[0]))
       self.dt=dt
       if len(self.initial_conditions)!=self.n:
           if len(self.initial_conditions)==1:
               self.initial_conditions=np.array([self.initial_conditions]*len(self.n))
           else:
               raise Exception("If the number 'n' of reservoirs is greater than 1 and initial conditions are not the same in these reservoirs, you must specify them. Check Boundaries")
            
    def executeRun(self):
        self.laggedInflow=np.zeros(len(self.Inflow)+self.lag)
        self.laggedInflow[self.lag:len(self.laggedInflow)]=self.Inflow
        self.routingSystem=LinearReservoirCascade(pars=[self.k,self.n],InitialConditions=self.initial_conditions,Boundaries=self.laggedInflow)
        self.routingSystem.computeOutFlow()
        self.Q=self.routingSystem.Outflow

#3. Modelos PQ/QQ

class HOSH4P1L(PydrologyProcedureInterface):
    """
    Modelo Operacional de Transformación de Precipitación en Escorrentía de 4 parámetros (estimables). Hidrología Operativa Síntesis de Hidrograma. Método NRCS, perfil de suelo con 2 reservorios de retención (sin efecto de base). Rutina de propagación por 'UH' arbitario (e.g. generado por triangularDistributon) o por función de transsferencia gamma. 
    """
    maxSurfaceStorage : float
    """Almacenamiento Máximo Superficial (reservorio de retención)"""
    maxSoilStorage : float
    """Almacenamiento máximo en el Suelo (reservorio de producción)"""
    Proc : str
    """Procedimiento de propagación ('Nash' o 'UH')"""
    Precipitation : np.ndarray
    """Precipitación (serie temporal)"""
    SurfaceStorage: np.ndarray
    """Almacenamiento en reservorio de retención (serie temporal)"""
    SoilStorage : np.ndarray
    """Almacenamiento en resservorio de producción (sserie temporal)"""
    Runoff : np.ndarray
    """Transferencia horizontal: escorrentía total (serie temporal)"""
    Infiltration : np.ndarray
    """Recarga de reservorio de producción (serie temporal)"""
    CumPrecip: np.ndarray
    "Precipitación acumulada durante el evento (serie temporal)"
    NetRainfall: np.ndarray
    """Precipitación neta (serie temporal)"""
    EVR1 : np.ndarray
    """Evapotranspiración real reservorio de abstracción"""
    EVR2 : np.ndarray
    """Evapotranspiración real reservorio de producción"""
    Q : float
    """Flujo encauzado (serie temporal)"""
    type='PQ Model'
    
    @property
    def maxSurfaceStorage(self) -> float:
        return self.pars[0]

    @property
    def maxSoilStorage(self) -> float:
        return self.pars[1]
    
    @property
    def u(self) -> float:
        if self.routingProc == 'UH':
            return self.pars[2]
        else:
            return None

    @property
    def k(self) -> float:
        if self.routingProc == 'Nash':
            return self.pars[2]
        else:
            return None
    
    @property
    def n(self) -> float:
        if self.routingProc == 'Nash':
            return self.pars[3]
        else:
            return None
    
    def __init__(self,pars : List[List[float]],Boundaries : Union[List[float],np.ndarray] =[[0],[0]],InitialConditions : Union[List[Tuple[float,float]],List[float]] =[0,0],Proc : str ='Nash'):
        """
            pars : List[List[float]]
                Lista con los valores de maxSurFaceStorage (reservorio de abstracción), maxSoilStorage (reservorio de producción) y parámetros tiempo de residencia (k) y n reservorios (caso Proc='Nash') o con último elemento como con ordenadas de Hidrograma Unitario (caso Proc='UH') 
            Boundaries : List[List[float]]
                Lista compuesta por listas de precipitación y evapotranspiración potencial (condiciones de borde) 
            InitialConditions : List[float]
                Lista con valores de almacenamiento inicial en reservorio de abstracción y en reservorio de producción
            Proc: str
                Procedimiento para transferencia: 'Nash' (cascada de Nash, debe proveerse k y n) o 'UH' (Hidrogramas Unitarios, array con j-vectores fila con valores de ordenadas)
        """
        super().__init__(pars,Boundaries,InitialConditions)
        self.routingProc=Proc
        self.soilSystem=SCSReservoirs(pars=[self.maxSurfaceStorage,self.maxSoilStorage])
        if self.routingProc == 'Nash':
                self.routingSystem=LinearChannel(pars=[self.k,self.n],Proc='Nash')    
        elif self.routingProc == 'UH':
                self.routingSystem=LinearChannel(pars=self.u,Proc='UH')
        else:
            raise Exception("invalid Proc. Must be one of: Nash, UH")
        self.Precipitation=np.array(self.boundaries[0],dtype='float')
        self.EVP=np.array(self.boundaries[1],dtype='float')
        self.EVR1=np.array([0]*len(self.Precipitation),dtype='float')
        self.EVR2=np.array([0]*len(self.Precipitation),dtype='float')
        self.SurfaceStorage=np.array([self.initial_conditions[0]]*(len(self.Precipitation)+1),dtype='float')
        self.SoilStorage=np.array([self.initial_conditions[1]]*(len(self.Precipitation)+1),dtype='float')
        self.NetRainfall=np.array([0]*len(self.Precipitation),dtype='float')
        self.Infiltration=np.array([0]*len(self.Precipitation),dtype='float')
        self.Runoff=np.array([0]*len(self.Precipitation),dtype='float')
        self.Q=np.array([0]*len(self.Precipitation),dtype='float')
    
    def computeRunoff(self): 
        j=0
        indexes=list()
        for row in list(self.Precipitation):
            if(self.Precipitation[j]>self.EVP[j]): #Condición de evento 
                indexes.append(j)
            else:
                if(len(indexes)>0): #Activa rutina de mojado (cómputo modelo de eventos SCS)
                    # logging.debug("ponding")
                    self.soilSystem.Precipitation=self.Precipitation[min(indexes):max(indexes)+1]-self.EVP[min(indexes):max(indexes)+1]
                    self.soilSystem.CumPrecip=np.array([0]*(len(self.soilSystem.Precipitation)),dtype='float')
                    self.soilSystem.NetRainfall=np.array([0]*(len(self.soilSystem.Precipitation)),dtype='float')
                    self.soilSystem.Infiltration=np.array([0]*(len(self.soilSystem.Precipitation)),dtype='float')
                    self.soilSystem.Runoff=np.array([0]*(len(self.soilSystem.Precipitation)),dtype='float')
                    self.soilSystem.SurfaceStorage=np.array([self.SurfaceStorage[min(indexes)]]*(len(self.soilSystem.Precipitation)+1),dtype='float')
                    self.soilSystem.SoilStorage=np.array([self.SoilStorage[min(indexes)]]*(len(self.soilSystem.Precipitation)+1),dtype='float')
                    self.soilSystem.computeAbstractionAndRunoff()
                    self.SurfaceStorage[min(indexes):max(indexes)+2]=self.soilSystem.SurfaceStorage
                    self.SoilStorage[min(indexes):max(indexes)+2]=self.soilSystem.SoilStorage
                    self.NetRainfall[min(indexes):max(indexes)+1]=self.soilSystem.NetRainfall
                    self.Infiltration[min(indexes):max(indexes)+1]=self.soilSystem.Infiltration
                    self.Runoff[min(indexes):max(indexes)+1]=self.soilSystem.Runoff
                    self.EVR1[min(indexes):max(indexes)+1]=self.EVP[min(indexes):max(indexes)+1]
                    indexes=list()
                if(len(indexes)==0): #Activa rutina de secado
                    # logging.debug("drying")
                    self.EVR1[j]=min(self.EVP[j],self.SurfaceStorage[j]+self.Precipitation[j])
                    self.NetRainfall[j]=max(0,self.Precipitation[j]-self.EVR1[j]+self.SurfaceStorage[j]-self.maxSurfaceStorage)
                    self.EVR2[j]=computeEVR(self.NetRainfall[j],self.EVP[j]-self.EVR1[j],self.SoilStorage[j],self.maxSoilStorage)
                    self.SurfaceStorage[j+1]=waterBalance(self.SurfaceStorage[j],self.Precipitation[j],self.EVR1[j]+self.NetRainfall[j])
                    self.Runoff[j]=max(0,self.NetRainfall[j]-self.EVR2[j]+self.SoilStorage[j-1]-self.maxSoilStorage)
                    self.SoilStorage[j+1]=waterBalance(self.SoilStorage[j],self.NetRainfall[j],self.EVR2[j]+self.Runoff[j])
            j=j+1            
   
    def computeOutFlow(self):
        self.routingSystem.Inflow=self.Runoff
        self.routingSystem.computeOutFlow()
        self.Q=self.routingSystem.Outflow
    
    def executeRun(self):
        self.computeRunoff()
        self.computeOutFlow()

class HOSH4P2L(PydrologyProcedureInterface):
    """
    Modelo Operacional de Transformación de Precipitación en Escorrentía de 4/6 parámetros (estimables), con 2 capas de suelo. Hidrología Operativa Síntesis de Hidrograma. Método NRCS, perfil de suelo con 2 reservorios de retención (zona superior) y un reservorio linear (zona inferior). Rutea utilizando una función respuesta de pulso unitario arbitraria o mediante na cascada de Nash (se debe especificar tiempo de residencia y número de reservorios)
    """
    maxSurfaceStorage : float
    """Almacenamiento Máximo Superficial (reservorio de retención)"""
    maxSoilStorage : float
    """Almacenamiento máximo en el Suelo (reservorio de producción)"""
    Proc : str
    """Procedimiento de propagación ('Nash' o 'UH')"""
    Precipitation : np.ndarray
    """Precipitación (serie temporal)"""
    SurfaceStorage: np.ndarray
    """Almacenamiento en reservorio de retención (serie temporal)"""
    SoilStorage : np.ndarray
    """Almacenamiento en resservorio de producción (sserie temporal)"""
    Runoff : np.ndarray
    """Transferencia horizontal: escorrentía total (serie temporal)"""
    Infiltration : np.ndarray
    """Recarga de reservorio de producción (serie temporal)"""
    CumPrecip: np.ndarray
    "Precipitación acumulada durante el evento (serie temporal)"
    NetRainfall: np.array
    """Precipitación neta (serie temporal)"""
    EVR1 : np.ndarray
    """Evapotranspiración real reservorio de abstracción"""
    EVR2 : np.ndarray
    """Evapotranspiración real reservorio de producción"""
    Q : float
    """Flujo encauzado (serie temporal)"""
    routingProc : str
    """Procedimiento de ruteo"""
    type='PQ Model'

    @property
    def maxSurfaceStorage(self) -> float:
        return self.pars[0]

    @property
    def maxSoilStorage(self) -> float:
        return self.pars[1]
    
    @property
    def phi(self) -> float:
        return self.pars[2]
    
    @property
    def kb(self) -> float:
        return self.pars[3]
    

    @property
    def tr(self) -> float:
        if self.routingProc == 'Nash':
            return self.pars[4]
        else:
            return None
        
    @property
    def n(self) -> float:
        if self.routingProc == 'Nash':
            return self.pars[5]
        else:
            return None
    
    @property
    def u(self) -> float:
        if self.routingProc == 'UH':
            return self.pars[4]
        else:
            return None

    def __init__(self,pars: List[List[float]],Boundaries : List[List[float]] =[[0],[0]],InitialConditions : List[float] =[0,0],Proc : str ='Nash'):
        """
            pars : List[List[float]]
                Lista con los valores de maxSurFaceStorage (reservorio de abstracción), maxSoilStorage (reservorio de producción), coeficiente de prorateo (flujo directo/flujo demorado, phi), coeficiente de recesión (autovalor, kb) y parámetros tiempo de residencia (k) y n reservorios (caso Proc='Nash') o con último elemento como lista incluyendo ordenadas de Hidrograma Unitario (caso Proc='UH') 
            Boundaries : List[List[float]]
                Lista o array2d compuesto por vectores (columna) de precipitación y evapotranspiración potencial (condiciones de borde) 
            InitialConditions : List[float]
                Lista con los valores de almacenamiento inicial en reservorio de abstracción y en reservorio de producción
            Proc: str
                Procedimiento para transferencia: 'Nash' (cascada de Nash, debe proveerse k y n) o 'UH' (Hidrogramas Unitarios, array con j-vectores fila con valores de ordenadas)
        """
        super().__init__(pars,Boundaries,InitialConditions)
        self.routingProc=Proc
        self.soilSystem=SCSReservoirs(pars=[self.maxSurfaceStorage,self.maxSoilStorage])
        if self.routingProc not in ['Nash','UH']:
            raise Exception("invalid Proc. Must be one of: Nash, UH")
        self.Precipitation=np.array(self.boundaries[0],dtype='float')
        self.EVP=np.array(self.boundaries[1],dtype='float')
        self.EVR1=np.array([0]*len(self.Precipitation),dtype='float')
        self.EVR2=np.array([0]*len(self.Precipitation),dtype='float')
        self.SurfaceStorage=np.array([self.initial_conditions[0]]*(len(self.Precipitation)+1),dtype='float')
        self.SoilStorage=np.array([self.initial_conditions[1]]*(len(self.Precipitation)+1),dtype='float')
        self.NetRainfall=np.array([0]*len(self.Precipitation),dtype='float')
        self.Infiltration=np.array([0]*len(self.Precipitation),dtype='float')
        self.Runoff=np.array([0]*len(self.Precipitation),dtype='float')
        self.Q=np.array([0]*len(self.Precipitation),dtype='float')
    def computeRunoff(self): 
        j=0
        indexes=list()
        for row in list(self.Precipitation):
            if(self.Precipitation[j]>self.EVP[j]): #Condición de evento 
                indexes.append(j)
            else:
                if(len(indexes)>0): #Activa rutina de mojado (cómputo modelo de eventos SCS)
                    # logging.debug("ponding")
                    self.soilSystem.Precipitation=self.Precipitation[min(indexes):max(indexes)+1]-self.EVP[min(indexes):max(indexes)+1]
                    self.soilSystem.CumPrecip=np.array([0]*(len(self.soilSystem.Precipitation)),dtype='float')
                    self.soilSystem.NetRainfall=np.array([0]*(len(self.soilSystem.Precipitation)),dtype='float')
                    self.soilSystem.Infiltration=np.array([0]*(len(self.soilSystem.Precipitation)),dtype='float')
                    self.soilSystem.Runoff=np.array([0]*(len(self.soilSystem.Precipitation)),dtype='float')
                    self.soilSystem.SurfaceStorage=np.array([self.SurfaceStorage[min(indexes)]]*(len(self.soilSystem.Precipitation)+1),dtype='float')
                    self.soilSystem.SoilStorage=np.array([self.SoilStorage[min(indexes)]]*(len(self.soilSystem.Precipitation)+1),dtype='float')
                    self.soilSystem.computeAbstractionAndRunoff()
                    self.SurfaceStorage[min(indexes):max(indexes)+2]=self.soilSystem.SurfaceStorage
                    self.SoilStorage[min(indexes):max(indexes)+2]=self.soilSystem.SoilStorage
                    self.NetRainfall[min(indexes):max(indexes)+1]=self.soilSystem.NetRainfall
                    self.Infiltration[min(indexes):max(indexes)+1]=self.soilSystem.Infiltration
                    self.Runoff[min(indexes):max(indexes)+1]=self.soilSystem.Runoff
                    self.EVR1[min(indexes):max(indexes)+1]=self.EVP[min(indexes):max(indexes)+1]
                    indexes=list()
                if(len(indexes)==0): #Activa rutina de secado
                    logging.debug("drying")
                    self.EVR1[j]=min(self.EVP[j],self.SurfaceStorage[j]+self.Precipitation[j])
                    self.NetRainfall[j]=max(0,self.Precipitation[j]-self.EVR1[j]+self.SurfaceStorage[j]-self.maxSurfaceStorage)
                    self.EVR2[j]=computeEVR(self.NetRainfall[j],self.EVP[j]-self.EVR1[j],self.SoilStorage[j],self.maxSoilStorage)
                    self.SurfaceStorage[j+1]=waterBalance(self.SurfaceStorage[j],self.Precipitation[j],self.EVR1[j]+self.NetRainfall[j])
                    self.Runoff[j]=max(0,self.NetRainfall[j]-self.EVR2[j]+self.SoilStorage[j-1]-self.maxSoilStorage)
                    self.SoilStorage[j+1]=waterBalance(self.SoilStorage[j],self.NetRainfall[j],self.EVR2[j]+self.Runoff[j])
            j=j+1                  
    def computeOutFlow(self):
        if self.routingProc == 'Nash':
            self.routingSystem=LinearChannel(pars=[self.tr,self.n],Boundaries=apportion(self.Runoff,self.phi))
        elif self.routingProc == 'UH':
            self.routingSystem=LinearChannel(pars=self.u,Boundaries=apportion(self.Runoff,self.phi),Proc='UH')
        else: 
            raise Exception("invalid Proc. Must be one of: Nash, UH")
        self.routingSystem.computeOutFlow()
        self.groundwaterSystem=LinearReservoirCascade(pars=[self.kb,1],Boundaries=apportion(self.Runoff,1-self.phi))
        self.groundwaterSystem.computeOutFlow()
        self.Q=self.routingSystem.Outflow[0:len(self.Runoff)]+self.groundwaterSystem.Outflow[0:len(self.Runoff)]
    def executeRun(self):
        self.computeRunoff()
        self.computeOutFlow()

class GR4J(PydrologyProcedureInterface):
    """
    Modelo Operacional de Transformación de Precipitación en Escorrentía de Ingeniería Rural de 4 parámetros (CEMAGREF). A diferencia de la versión original, la convolución se realiza mediante producto de matrices. Parámetros: Máximo almacenamiento en reservorio de producción, tiempo al pico (hidrograma unitario),máximo almacenamiento en reservorio de propagación, coeficiente de intercambio.
    """
    Precipitation : np.ndarray
    """Precipitación (serie temporal)"""
    EVP : np.ndarray
    """Evapotranspiración (serie temporal)"""
    Runoff : np.ndarray
    """Escorentía reservorio de producción (serie temporal)"""
    Q: np.ndarray
    """Flujo encauzado (serie temporal)"""
    RoutingProc : str
    """Procedimiento de ruteo"""
    type='PQ Model'

    @property
    def prodStoreMaxStorage(self) -> float:
        return self.pars[0]

    @property
    def T(self) -> float:
        return self.pars[1]
    
    @property
    def routStoreMaxStorage(self) -> float:
        return self.pars[2]

    @property
    def waterExchange(self) -> float:
        if len(self.pars)<4:
            return 0
        else:
            return self.pars[3]    

    def __init__(self,pars : List[float], Boundaries : List[List[float]] =[[0],[0]],InitialConditions : List[float]=[0,0],Proc='CEMAGREF SH'):
        """
            pars : List[float]
                Lista con los valores de Almacenamiento Máximo en Reservorio de Producción, Tiempo al pico , Almacenamiento Máximo en Reservorio de Tránsito y coeficiente de intercambio o fugas 
            Boundaries : List[List[float]]
                Lista con series de precipitación y  de evapotranspiración potencial, declaradas como litas (pmad : float, etpd : float) (condiciones de borde) 
            InitialConditions : List[float]
                Lista con valores de almacenamiento inicial en reservorio de producción y en reservorio de tránsito
        """
        super().__init__(pars,Boundaries,InitialConditions)
        self.RoutingProc=Proc
        self.u1=grXDistribution(self.T,distribution='SH1')
        self.u2=grXDistribution(self.T,distribution='SH2')
        self.Precipitation=np.array(self.boundaries[0],dtype='float')
        self.EVP=np.array(self.boundaries[1],dtype='float')
        self.Runoff=np.array([0]*len(self.Precipitation),dtype='float')
        self.Q=np.array([0]*len(self.Precipitation),dtype='float')
        self.prodStore=ProductionStoreGR4J(pars=[self.prodStoreMaxStorage],Boundaries=[self.Precipitation,self.EVP],InitialConditions=[self.initial_conditions[0]]) #cambios realizados en ProdStore (Ajuste interfaz) 20240524
    def computeRunoff(self):
        self.prodStore.computeOutFlow()
        self.Runoff=self.prodStore.Runoff
    def computeOutFlow(self):
        self.channel1=LinearChannel(pars=self.u1,Boundaries=apportion(0.9,self.Runoff),Proc='UH')
        self.channel2=LinearChannel(pars=self.u2,Boundaries=apportion(0.1,self.Runoff),Proc='UH')
        self.channel1.computeOutFlow()
        self.channel2.computeOutFlow()
        self.routStore=RoutingStoreGR4J(pars=[self.routStoreMaxStorage,self.waterExchange],Boundaries=[self.channel1.Outflow],InitialConditions=[self.initial_conditions[1]]) #cambios realizados en RoutStore 20240524. Hay que ajustar todos los procedimientos a los requisitos impuestos por interfaz, continuar
        self.routStore.computeOutFlow()
        n=min(len(self.routStore.Runoff),len(self.channel2.Outflow))
        self.DirectRunoff=np.array([0]*n,dtype='float')
        for j in range(0,n):
            self.DirectRunoff[j]=max(0,self.channel2.Outflow[j]+self.routStore.Leakages[j])
        self.Q=self.routStore.Runoff[0:n]+self.DirectRunoff[0:n]
    def executeRun(self):
        self.computeRunoff()
        self.computeOutFlow()

class HIDROSAT(PydrologyProcedureInterface):
    """
    Modelo Operacional de Transformación de Precipitación en Escorrentía HIDROSAT. 
    """
    S0: float 
    """Capacidad de campo"""
    W0 : float
    """Almacenamiento de referencia en planicie aluvial"""
    Q0 : float
    """Caudal de referencia"""
    gamma: float
    """Factor de forma (abstracción) en relación Q(W)"""
    epsilon: float
    """Error máximo tolerable en Newthon-Raphson"""
    dt : float
    """Longitud del paso de cálculo"""
    Precipitation : np.ndarray
    """Precipitación Media Areal"""
    EVP: np.ndarray
    """Evapotranspiración Potencial Media Areal"""
    soilStorage: np.ndarray
    """Almacenamiento en reservorio de retención (agua de tensión/interfluvios)"""
    EVSoil: np.ndarray
    """Evapotranspiración (retención interfluvios)"""
    freeWater: np.ndarray
    """Agua gravífica"""
    DirectRunoff: np.ndarray
    """Escorrentía directa"""
    Runoff: np.ndarray
    """Escorrentía demorada"""
    floodplainStorage: np.ndarray
    """Almacenamiento en reservorio de detención (propagación escorrentía/planicie aluvial)"""
    EVFloodPlain: np.ndarray
    """Evaporación (detención planicie aluvial)"""
    Q: np.ndarray
    "Caudal"
    type='HIDROSAT Model (Giordano, 2014)'

    @property
    def S0(self) -> float:
        return self.pars[0]
    
    @property
    def K(self) -> float:
        return self.pars[1]
    
    @property
    def N(self) -> float:
        return self.pars[2]
    
    @property
    def W0(self) -> float:
        return self.pars[3]
    
    @property
    def Q0(self) -> float:
        return self.pars[4]
    
    @property
    def gamma(self) -> float:
        return self.pars[5]

    @property
    def maxFlooded(self) -> float:
        if(len(self.pars)<7):
            return 1
        else:
            if self.pars[6]>1 or self.pars[6]<0:
                raise NameError("Max Flooded Area is a fraction of basin area, between 0-1")
            else:
                return self.pars[6]
    
    @property
    def detentionRatio(self) -> float:
        if(len(self.pars)<8):
            return 1
        else:
            if self.pars[7]>1 or self.pars[7]<0:
                raise NameError("Max Flooded Area is a fraction of basin area, between 0-1")
            else:
                return self.pars[7]

    @property
    def epsilon(self) -> float:
        if(len(self.pars)<9):
            return 0.0001
        else:
            return self.pars[8]
    
    def __init__(self,pars : List[float], Boundaries : List[List[float]] =[[0],[0],[0]],InitialConditions : List[float]=[0,0,0,0],dt : float = 1):
        """
            pars : List[float]
                Lista con los valores de capacidad de campo (almacenamiento máximo reservorio de retención), tiempo de residencia, cantidad de reservorios lineales (flujo demorado), almacenamiento de referencia en planicie aluvial (almacenamiento máximo reservorio de detención), caudal de referencia (caudal máximo reservorio de detención), factor de forma (fator de escala entre AS/AS0 y Q/Q), máxima área anegada (si no asume 1) y opcionalmente detentionRatio (proporción de flujo directo interceptado por reservorio de detención) y epsilon (parámetro de corte en Newton Raphson) 
            Boundaries : List[List[float]]
                Lista con series de precipitación y  de evapotranspiración potencial, y opcionalmente de aporte desde aguas arriba, declaradas como listas (pmad : float, etpd : float, inflow : float) (condiciones de borde) 
            InitialConditions : List[float]
                Lista con valores de almacenamiento inicial en perfil de suelo, caudal inicial en cascada de reservorios lineales, almacenamiento y área anegada iniciales en planicie aluvial
            dt : float
                Longitud del subpaso de cómputo en Newton-Raphson  
        """
        super().__init__(pars,Boundaries,InitialConditions)
        self.Precipitation=np.array(self.boundaries[0],dtype='float')
        self.EVP=np.array(self.boundaries[1],dtype='float')
        if len(self.boundaries) < 3:
            self.inFlow=np.array([0]*len(self.Precipitation),dtype='float')
        else:
            self.inFlow=np.array(self.boundaries[2],dtype='float')
        self.soilStorage=np.array([self.initial_conditions[0]]*(len(self.Precipitation)+1),dtype='float')
        self.EVSoil=np.array([0]*len(self.Precipitation),dtype='float')
        self.freeWater=np.array([0]*len(self.Precipitation),dtype='float')
        self.Runoff=np.array([self.initial_conditions[1]]*(len(self.Precipitation)+1),dtype='float')
        self.DirectRunoff=np.array([0]*len(self.Precipitation),dtype='float')
        self.floodplainStorage=np.array([self.initial_conditions[2]]*(len(self.Precipitation)+1),dtype='float')
        self.EVFloodPlain=np.array([0]*len(self.Precipitation),dtype='float')
        self.Q=np.array([self.Q0*(self.initial_conditions[2]/self.W0)**(self.gamma)]*(len(self.Precipitation)+1),dtype='float')
        self.dt=dt
        if len(self.initial_conditions) < 4:
            self.Flooded=np.array([0]*(len(self.Precipitation)+1),dtype='float')
        else:
            self.Flooded=np.array([self.initial_conditions[3]]*(len(self.Precipitation)+1),dtype='float')
        
    def executeRun(self):
        for i in range(0,len(self.Precipitation)-1,1):
            phi=(1-self.Flooded[i])
            p_soil=[phi*self.Precipitation[i]]
            evp_soil=[phi*self.EVP[i]]
            initial_soil=[self.soilStorage[i]]
            self.soilSystem=RetentionReservoir(pars=[self.S0],Boundaries=[p_soil,evp_soil],InitialConditions=initial_soil)
            self.soilSystem.computeRunoff()
            self.EVSoil[i]=self.soilSystem.EV[0]
            self.freeWater[i]=self.soilSystem.Runoff[0]
            self.soilStorage[i+1]=self.soilSystem.Storage[1]
            inflow_gw=[self.freeWater[i]]
            if i==0:
                initial_gw=[self.initial_conditions[1]]
            else:
                initial_gw=[self.gwSystem.Cascade[0,0],self.gwSystem.Cascade[0,1],self.gwSystem.Cascade[1,0],self.gwSystem.Cascade[1,1]]
            self.gwSystem=LinearReservoirCascade(pars=[self.K,self.N],Boundaries=inflow_gw,InitialConditions=initial_gw)
            self.gwSystem.computeOutFlow()
            self.Runoff[i+1]=self.gwSystem.Outflow[1]
            self.DirectRunoff[i]=(1-phi)*self.Precipitation[i]
            self.EVFloodPlain[i]=min(self.floodplainStorage[i]+self.DirectRunoff[i],(1-phi)*self.EVP[i])
            p_floodplain=[apportion(self.DirectRunoff[i],self.detentionRatio)]
            q_direct=[apportion(self.DirectRunoff[i],1-self.detentionRatio)]
            ev_floodplain=[self.EVFloodPlain[i]]
            inflows_floodplain=[self.Runoff[i]+self.inFlow[i],self.Runoff[i+1]+self.inFlow[i+1]]
            initial_floodplain=[self.floodplainStorage[i],self.Flooded[i]]
            self.routingSystem=HIDROSATPowerLawReservoir(pars=[self.W0,self.Q0,self.gamma,self.maxFlooded,self.epsilon],Boundaries=[p_floodplain,ev_floodplain,inflows_floodplain[0],inflows_floodplain[1]],InitialConditions=initial_floodplain,dt=self.dt)
            self.routingSystem.computeOutFlow()
            self.floodplainStorage[i+1]=self.routingSystem.Storage[1]
            self.Flooded[i+1]=self.routingSystem.Flooded[1]
            self.Q[i+1]=self.routingSystem.Outflow[1]+q_direct

if __name__ == "__main__":
    import sys


    