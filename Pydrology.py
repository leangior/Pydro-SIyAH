#/usr/bin/python3
#Librería de métodos para modelación hidrológica SSIyAH-INA, 2022
import math
from sys import maxsize
from zlib import MAX_WBITS
import numpy  as np
import pandas as pd
import matplotlib.pyplot as plt
import os, glob

#genera gráfica de señales de entrada/salida para diagnóstico visual

def testPlot(Inflow,Outflow):
    plt.plot(Inflow,'r')
    plt.plot(Outflow,'b')
    plt.show()

#desplaza subíndice de serie una unidad a la izquierda (ajuste para representación discreta)

def shiftLeft(array1d,fill=0):
    i=0
    shift=np.array([0]*len(array1d),dtype='float')
    for value in range(1,len(array1d)):
        shift[i]=array1d[i+1]
        i=i+1
    shift[i]=fill
    return(shift)

#Importa CSV a estructura de datos de entrada PQ diario (t,PMA,EV0)

def getDrivers(file,tCol='t'):
    data=pd.read_csv(file)
    data[tCol]=pd.to_datetime(data[tCol],format='%Y-%m-%d')
    data.index=data[tCol]
    return(data)

#Crea Condiciones de Borde p y evp para PQ 
def makeBoundaries(p=[0],evp=[0]):
    boundaries=np.array([[0]*2]*len(p),dtype='float')
    boundaries[:,0]=p
    boundaries[:,1]=evp
    return boundaries

#Diferencia una serie
def differentiate(list,asume_initial='yes'):
    if(asume_initial=='yes'):
        dif=[list[0]]*len(list)
    else:
        dif=[0]*len(list)
    for i in range(1,len(list)):
        dif[i]=list[i]-list[i-1]
    return dif

#Integra por método del trapecio
def integrate(list,dt):
    int=0
    for i in range(1,len(list)):
        int=int+(list[i]+list[i-1])*dt/2
    return int

#Computa Hidrogramas Triangulares (Función Respuesta Unitaria, método Triangular Simétrico o Triangular SCS)
def triangularDistribution(T,distribution='Symmetric',dt=0.01,shift='T',approx='T'):
    if distribution == 'Symmetric':
        tb=2*T
        peakValue=1/T
    if distribution == 'SCS':
        tb=8/3*T
        peakValue=3/4*1/T
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
    if approx=='T':
        U=U/sum(U)
    if shift == 'T':
        U=shiftLeft(U)         
        return(U[0:(len(U)-1)])
    else:
        return(U)

#Computa Función Respuesta Unitaria Cascada de n reservorios Lineales con tiempo de residencia k, obtenida por integración numérica a resolución dt (método del trapecio). El parámetro shift se agregó para desplazar los subíndices una unidad a la izquierda, puesto que si no muestrea la integración a fin de intervalo de cómputo, pudiéndose introducir un artefacto numérico con efecto de retardo, aproximadamente en una unidad.
def gammaDistribution(n,k,dt=1,m=10,approx='T',shift='T'):
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
    if approx == 'T':
        U=U/sum(U)
    if shift == 'T':
        U=shiftLeft(U)
    return U

#Computa HUs propuestos en modelos GRX (GR4J, GRP) 
def grXDistribution(T,distribution='SH1',dt=0.5,approx='T',Agg='T'):    
    if distribution == 'SH1':
        tb=T
        k=1
    if distribution == 'SH2':
        tb=2*T
        k=1/2
    ndimu=int(round(tb/dt,0)+2)
    ndimU=int(round(tb,0)+1)
    u=np.array([0]*(ndimu),dtype='float')
    U=np.array([0]*(ndimU),dtype='float')
    for t in np.array(list(range(0,ndimu))):
        if t*dt<T:
            u[t]=k*(t*dt/T)**(5/2)
        if t*dt>T and t*dt<tb:
             u[t]=1-k*(2-t*dt/T)**(5/2)        
        else:
           if t*dt>tb:
                u[t]=1
    u=differentiate(u,'no')
    for j in range(0,ndimU):
        min=int(j/dt)
        max=int((j+1)/dt)
        U[j]=sum(u[min:max])
    if Agg == 'T':
        return(U)
    else:
        return(u)
        
#Computa Matriz de pulsos para Convolución con At 12:00 on day-of-month 1.”
def getPulseMatrix(inflows,u):
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

#Computa Ecuación de Conservación
def waterBalance(Storage=0,Inflow=0,Outflow=0):
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

#Realiza correción de sesgo por simple updating (método propuesto por el Servicio Ruso)
def SimonovKhristoforov(sim,obs): 
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
class RetentionReservoir:
    """
    Reservorio de Retención. Vector pars de sólo parámetro: capacidad máxima de abstracción [MaxStorage]. Condiciones Iniciales (InitialConditions): [Initial Storage] (ingresa como vector). Condiciones de Borde (Boundaries): vectior [[Inflow],[EV]]. 
    """
    type='Retention Reservoir'
    def __init__(self,pars,InitialConditions=[0],Boundaries=[[0],[0]],Proc='Abstraction'):
        self.MaxStorage=pars[0]
        self.Inflow=np.array(Boundaries[0],dtype='float')
        self.EV=np.array(Boundaries[1],dtype='float')
        self.Storage=np.array([InitialConditions[0]]*(len(self.Inflow)+1),dtype='float')
        self.Runoff=np.array([0]*len(self.Inflow),dtype='float')
        self.Proc=Proc
        self.dt=1
    def computeRunoff(self):
        for i in range(0,len(self.Inflow),1):
            if self.Proc == 'Abstraction':
                self.Runoff[i]=max(0,self.Inflow[i]-self.EV[i]+self.Storage[i]-self.MaxStorage)
            if self.Proc == 'CN_h0_continuous':
                self.Runoff[i]=(max(self.Inflow[i]-self.EV[i],0))**2/(self.MaxStorage-self.Storage[i]+self.Inflow[i]-self.EV[i])
            self.Storage[i+1]=waterBalance(self.Storage[i],self.Inflow[i],self.EV[i]+self.Runoff[i])


#1.B Reservorio Lineal. ###REVISAR RUTINA TOMA LEN() OF UNSIZED OBJECT 
class LinearReservoir:
    """
    Reservorio Lineal. Vector pars de un sólo parámetro: Tiempo de residencia (K). Vector de Condiciones Iniciales (InitialConditions): Storage, con el cual computa Outflow. Condiciones de Borde (Boundaries): Inflow y EV.
    """
    type='Linear Reservoir'
    def __init__(self,pars,InitialConditions=[0],Boundaries=[[0],[0]],Proc='Agg',dt=1):
        self.K=pars[0]
        self.Inflow=np.array(Boundaries[0],dtype='float') 
        self.EV=np.array(Boundaries[1],dtype='float')
        self.Storage=np.array([InitialConditions[0]]*(len(self.Inflow)+1),dtype='float')
        self.Outflow=(1/self.K)*self.Storage
        self.Proc=Proc
        if Proc == ('Agg' or 'API'):
            self.dt=1
        if Proc == 'Instant':
            self.dt=dt
    def computeOutFlow(self):
        for i in range (0,len(self.Inflow),1):
            if self.Proc == 'Agg':
                self.Outflow[i]=(1/self.K)*(self.Storage[i]+self.Inflow[i])
                self.Storage[i+1]=waterBalance(self.Storage[i],self.Inflow[i],self.Outflow[i])
            if self.Proc == 'Instant':
                end=int(1/self.dt+1)
                Storage=self.Storage[i]
                Outflow=self.Outflow[i]    
                for t in range(1,end,1):
                    Storage=waterBalance(Storage,self.Inflow[i]*self.dt,Outflow*self.dt)
                    Outflow=(1/self.K)*(Storage)
                self.Storage[i+1]=Storage
                self.Outflow[i+1]=Outflow
            if self.Proc == 'API':
                self.Storage[i+1]=(1-1/self.K)*self.Storage[i]+self.Inflow[i]
                self.Outflow[i+1]=(1/self.K)*self.Storage
                
class ProductionStoreGR4J:
    """
    Reservorio de Producción de Escorrentía modelo GR4J
    """
    type='GR4J Runoff Production Store'
    def __init__(self,pars,InitialConditions=[0],Boundaries=[[0],[0]],Proc='Time Discrete Agg'):
        self.MaxSoilStorage=pars[0]
        self.Precipitation=np.array(Boundaries[:,0],dtype='float')
        self.EVP=np.array(Boundaries[:,1],dtype='float')
        self.SoilStorage=np.array([InitialConditions[0]]*(len(self.Precipitation)+1),dtype='float')
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
            self.Recharge[i]=(self.MaxSoilStorage*(1-(relativeMoisture)**2)*np.tanh(ratio_netRainfall_maxStorage))/(1+relativeMoisture*ratio_netRainfall_maxStorage)
            self.EVR[i]=(self.SoilStorage[i]*(2-relativeMoisture)*np.tanh(ratio_netEVP_maxStorage))/(1+(1-relativeMoisture)*np.tanh(ratio_netEVP_maxStorage))
            self.SoilStorage[i+1]=waterBalance(self.SoilStorage[i],self.Recharge[i],self.EVR[i])
            self.Infiltration[i]=self.MaxSoilStorage*(1-((1+4/9*relativeMoisture)**4)**(-1/4))
            self.SoilStorage[i+1]=waterBalance(self.SoilStorage[i+1],0,self.Infiltration[i])
            self.Runoff[i]=self.Infiltration[i]+self.NetRainfall[i]-self.Recharge[i]

class RoutingStoreGR4J:
    """
    Reservorio de Propagación de Escorrentía modelo GR4J
    """
    type='GR4J Runoff Routing Store'
    def __init__(self,pars,InitialConditions=[0],Boundaries=[0],Proc='Time Discrete Agg'):
        self.MaxStorage=pars[0]
        if not pars[1]:
            self.waterExchange=0
        else:
            self.waterExchange=pars[1]
        self.Inflow=np.array(Boundaries,dtype='float')
        self.Leakages=np.array([0]*len(self.Inflow),dtype='float')
        self.Runoff=np.array([0]*len(self.Inflow),dtype='float')
        self.Storage=np.array([InitialConditions[0]]*(len(self.Inflow)+1),dtype='float')
    def computeOutFlow(self):
         for i in range(0,len(self.Inflow)):
            relativeMoisture=self.Storage[i]/self.MaxStorage
            self.Leakages[i]=self.waterExchange*relativeMoisture**(7/2)
            self.Storage[i+1]=max(0,self.Storage[i]+self.Inflow[i])
            relativeMoisture=self.Storage[i+1]/self.MaxStorage
            self.Runoff[i]=self.Storage[i+1]*(1-(1+relativeMoisture**4)**(-1/4))
            self.Storage[i+1]=waterBalance(self.Storage[i+1],0,self.Runoff[i])

class SCSReservoirs:
    """
    Sistema de 2 reservorios de retención (intercepción/abstracción superficial y retención en perfil de suelo - i.e. capacidad de campo-), con función de cómputo de escorrentía siguiendo el método propuesto por el Soil Conservation Service. Vector pars de dos parámetros: Máximo Almacenamiento Superficial (Abstraction) y Máximo Almacenamiento por Retención en Perfil de Suelo (MaxStorage). Condiciones iniciales: Almacenamiento Superficial y Almacenamiento en Perfil de Suelo (lista de valores). Condiciones de Borde: Hietograma (lista de valores).
    """
    type='Soil Conservation Service Model for Runoff Computation (Curve Number Method / Discrete Approach)'
    def __init__(self,pars,InitialConditions=[0,0],Boundaries=[0],Proc='Time Discrete Agg'):
        self.MaxSurfaceStorage=pars[0]
        self.MaxStorage=pars[1]
        self.Precipitation=np.array(Boundaries,dtype='float')
        self.SurfaceStorage=np.array([InitialConditions[0]]*(len(self.Precipitation)+1),dtype='float')
        self.SoilStorage=np.array([InitialConditions[1]]*(len(self.Precipitation)+1),dtype='float')
        self.Runoff=np.array([0]*len(self.Precipitation),dtype='float')
        self.Infiltration=np.array([0]*len(self.Precipitation),dtype='float')
        self.CumPrecip=np.array([0]*len(self.Precipitation),dtype='float')
        self.NetRainfall=np.array([0]*len(self.Precipitation),dtype='float') 
        self.Proc=Proc
        self.dt=1
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
            self.SurfaceStorage[i+1]=min(self.SurfaceStorage[0]+Abstraction,self.CumPrecip[i])
        self.Runoff=differentiate(self.Runoff)
        self.NetRainfall=differentiate(self.NetRainfall)
        self.Infiltration=differentiate(self.Infiltration)
        for i in range(0,len(self.SoilStorage)-1):
            self.SoilStorage[i+1]=waterBalance(self.SoilStorage[i],self.Infiltration[i])        

class SCSReservoirsMod:
    """
    Sistema de 2 reservorios de retención+detención (una capa de abstracción superficial/suelo y otra capa de retención/detención en resto perfil de suelo), con función de cómputo de escorrentía siguiendo el método propuesto por el Soil Conservation Service y añadiendo pérdida continua por flujo de base (primario). Vector pars de 3 parámetros: Máxima Abtracción por retención (Abstraction) y Máximo Almacenamiento por Retención+Detención en Perfil de Suelo (MaxStorage) y coefiente de pérdida K. Se añade pérdida continua. Condiciones iniciales: Almacenamiento Superficial y Almacenamiento en Perfil de Suelo (lista de valores). Condiciones de Borde: Hietograma (lista de valores).
    """
    type='Soil Conservation Service Model for Runoff Computation (Curve Number Method / Discrete Approach)'
    def __init__(self,pars,InitialConditions=[0,0],Boundaries=[0],Proc='Time Discrete Agg'):
        self.MaxSurfaceStorage=pars[0]
        self.MaxStorage=pars[1]
        self.K=pars[2]
        self.Precipitation=np.array(Boundaries,dtype='float')
        self.SurfaceStorage=np.array([InitialConditions[0]]*(len(self.Precipitation)+1),dtype='float')
        self.SoilStorage=np.array([InitialConditions[1]]*(len(self.Precipitation)+1),dtype='float')
        self.Runoff=np.array([0]*len(self.Precipitation),dtype='float')
        self.Infiltration=np.array([0]*len(self.Precipitation),dtype='float')
        self.CumPrecip=np.array([0]*len(self.Precipitation),dtype='float')
        self.NetRainfall=np.array([0]*len(self.Precipitation),dtype='float')
        self.BaseFlow=np.array([0]*len(self.Precipitation),dtype='float') 
        self.Proc=Proc
        self.dt=1
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
            self.SurfaceStorage[i+1]=min(self.SurfaceStorage[0]+Abstraction,self.CumPrecip[i])
        self.Runoff=differentiate(self.Runoff)
        self.NetRainfall=differentiate(self.NetRainfall)
        self.Infiltration=differentiate(self.Infiltration)
        for i in range(0,len(self.SoilStorage)-1):
            self.BaseFlow[i]=(1-self.K)*(self.SoilStorage[i]+self.Infiltration[i])
            self.SoilStorage[i+1]=waterBalance(self.SoilStorage[i],self.Infiltration[i],self.BaseFlow[i])        



#2. Proceso Q-Q: Componentes de Función Distribución de Escorrentía o Tránsito Hidrológico

#2.A Cascada de Reservorios Lineales (Discreta). Dos parámetros: Tiempo de Resdiencia (K) y Número de Reservorios (N)
class LinearReservoirCascade:
    """
    Cascada de Reservorios Lineales (Discreta). Vector pars de dos parámetros: Tiempo de Residencia (K) y Número de Reservorios (N). Vector de Condiciones Iniciales (InitialConditions): Si es un escalar (debe ingresarse como elemento de lista) genera una matriz de 2xN con valor constante igual al escalar, también puede ingresarse una matriz de 2XN que represente el caudal inicial en cada reservorio de la cascada. Condiciones de Borde (Boundaries): vector Inflow. 
    """
    type='Discrete Cascade of N Linear Reservoirs with Time K'
    def __init__(self,pars,Boundaries=[0],InitialConditions=[0],create='yes',Proc='Discretely Coincident',dt=1):
        self.K=pars[0]
        if not pars[1]:
            self.N=2
        else:
            self.N=pars[1]
        self.Inflow=np.array(Boundaries)   
        if  create == 'yes':
            self.Cascade=np.array([[InitialConditions[0]]*self.N]*2,dtype='float')
        else:
            self.Cascade=np.array(InitialConditions,dtype='float')
        self.Outflow=np.array([InitialConditions[0]]*(len(Boundaries)+1),dtype='float')
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

#2.B Canal Muskingum 
# EN DESARROLLO (MUSKINGUM y CUNGE) --> VER RESTRICCIONES NUMÉRICAS y SI CONSIDERAR CURVAS HQ y BH COMO PARAMETROS DEL METODO. POR AHORA FINALIZADO MUSKINGUM CLÁSICO. CUNGE DEBE APOYARSE SOBRE EL MISMO, MODIFICANDO PARS K y X
class MuskingumChannel:
    """
    Método de tránsito hidrológico de la Oficina del río Muskingum. Vector pars de dos parámetros: Tiempo de Tránsito (K) y Factor de forma (X) [Proc='Muskingum'] o . Condiciones Iniciales (InitialConditions): matriz de condiciones iniciales o valor escalar constante. Condiciones de borde: Hidrograma en nodo superior de tramo. 
    """
    #A fin de mantener condiciones de estabilidad numérica en la propagación (conservar volumen), sobre la base de la restricción 2KX<=dt<=2K(1-X) (Chin,2000) y como dt viene fijo por la condición de borde (e.g. por defecto 'una unidad') y además se pretende respetar el valor de K, se propone incrementar la resolución espacial dividiendo el tramo en N subtramos de igual longitud, con tiempo de residencia mínimo T=K/N, para el caso dt<2KX (frecuencia de muestreo demasiado alta). Luego, aplicando el criterio de chin se sabe que el valor crítico de dt debe satisfacer dt=uT, específicamente con u=2X y T = K/N--> N=2KX/dt. Al mismo tiempo si dt>2K(1-X) (frecuencia de muestreo demasiado baja), el paso de cálculo se subdivide en M subpasos de longitud dT=2K(1-X) de forma tal que dT/dt=dv y M=dt/dv. Self.tau especifica el subpaso de cálculo (siendo self.M la cantidad de subintervalos utilizados) y self.N la cantidad de subtramos. 
    type='Muskingum Channel'
    def __init__(self,pars,Boundaries=[0],InitialConditions=[0],Proc='Muskingum Routing Method',dt=1):
        self.K=pars[0]
        self.X=pars[1]
        self.dt=dt
        self.lowerbound=2*self.K*self.X
        self.upperbound=2*self.K*(1-self.X)
        self.Inflow=np.array(Boundaries,dtype='float')
        self.Outflow=np.array([0]*len(self.Inflow),dtype='float')
        self.InitialConditions=np.array(InitialConditions,dtype='float')
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
            

#2.C Tránsito Lineal con funciones de transferencia. Por defecto, se asume una distrinución gamma con parámetros n (número de reservorios) y k (tiempo de residencia). Asimismo, se considera n=2, de modo tal que tp=k (el tiempo al pico es igual al tiempo de residencia) 
class LinearChannel:
    """
    Método de tránsito hidrológico implementado sobre la base de teoría de sistemas lineales. Así, considera al tránsito de energía, materia o información como un proceso lineal desde un nodo superior hacia un nodo inferior. Específicamente, sea I=[I1,I2,...,IN] el vector de pulsos generados por el borde superior y U=[U1,U2,..,UM] una función de distribución que representa el prorateo de un pulso unitario durante el tránsito desde un nodo superior (borde) hacia un nodo inferior (salida), el sistema opera aplicando las propiedades de proporcionalidad y aditividad, de manera tal que es posible propagar cada pulso a partir de U y luego mediante la suma de estos prorateos obtener el aporte de este tránsito sobre el nodo inferior (convolución).
    """
    type='Single Linear Channel'
    def __init__(self,pars,Boundaries=[0],Proc='Nash',dt=1):
       self.pars=np.array(pars,dtype='float')
       self.Inflow=np.array(Boundaries,dtype='float')
       self.Proc=Proc
       self.dt=dt
       if self.Proc == 'Nash':
            self.k=self.pars[0]
            self.n=self.pars[1]
            self.u=gammaDistribution(self.n,self.k,self.dt)
       if self.Proc == 'UH':
            self.u=self.pars
       self.Outflow=np.array([[0]]*(len(self.Inflow)+len(self.u)-1))
    def computeOutFlow(self):
        I=getPulseMatrix(self.Inflow,self.u)
        self.Outflow=np.dot(I,self.u)

class LinearNet:
    """
    Método de tránsito hidrológico implementado sobre la base de teoría de sistemas lineales. Así, considera al tránsito de energía, materia o información como un proceso lineal desde N nodos superiores hacia un nodo inferior. Específicamente, sea I=[I1,I2,...,IN] un vector de pulsos generados por un borde y U=[U1,U2,..,UM] una función de distribución que representa el prorateo de un pulso unitario durante el tránsito desde un nodo superior (borde) hacia un nodo inferior (salida), aplicando las propiedades de proporcionalidad y aditividad es posible propagar cada pulso a partir de U y luego mediante su suma obtener el aporte de este tránsito sobre el nodo inferior, mediante convolución. Numéricamente el sistema se representa como una transformación matricial (matriz de pulsos*u=vector de aportes). Consecuentemente, el tránsito se realiza para cada borde y la suma total de estos tránsitos constituye la señal transitada sobre el nodo inferior.  Condiciones de borde: array 2D con hidrogramas en nodos superiores del tramo, por columna. Parámetros: función de distribución (proc='EmpDist') o tiempo de residencia (k) y número de reservorios (n), si se desea utilizar el método de hidrograma unitario de Nash (proc='Nash'), pars es un array bidimensional en donde la información necesaria para cada nodo se presenta por fila (parámetros de nodo). El parámetro dt refiere a la longitud de paso de cálculo para el método de integración, siendo dt=1 la resolución nativa de los hidrogramas de entrada provistos. Importante, las funciones de transferencia deben tener la misma cantidad de ordenadas (dimensión del vector) 
    """
    type='Linear Routing System. System of Linear Channels'
    def __init__(self,pars,Boundaries=[0],Proc='Nash',dt=1):
        self.pars=np.array(pars,dtype='float')
        self.Inflows=np.array(Boundaries,dtype='float')
        self.Proc=Proc
        self.dt=dt
    def computeOutflow(self):
        j=0
        for channel_j in range(1,len(self.Inflows[0,:])+1):
            linear=LinearChannel(pars=self.pars[j,:],Boundaries=self.Inflows[:,j],dt=self.dt)
            linear.computeOutFlow()
            if j==0:
                self.Outflow=linear.Outflow
            if j>0:
                nrows=max(len(self.Outflow),len(linear.Outflow))
                f=np.zeros((nrows))
                if len(self.Outflow) > len(linear.Outflow):
                   f[0:len(linear.Outflow)]=linear.Outflow[0:len(linear.Outflow)]
                   self.Outflow=self.Outflow+f 
                else:
                   f[0:len(self.Outflow)]=self.Outflow[0:len(self.Outflow)]
                   self.Outflow=f+linear.Outflow 
            j=j+1

#3. Modelos PQ/QQ

class HOSH4P1L:
    """
    Modelo Operacional de Transformación de Precipitación en Escorrentía de 4 parámetros (estimables). Hidrología Operativa Síntesis de Hidrograma. Método NRCS, perfil de suelo con 2 reservorios de retención (sin efecto de base).
    """
    type='PQ Model'
    def __init__(self,pars,Boundaries=[0],InitialConditions=[[0],[0]],Proc='Nash'):
        self.maxSurfaceStorage=pars[0]
        self.maxSoilStorage=pars[1]
        self.soilSystem=SCSReservoirs(pars=[self.maxSurfaceStorage,self.maxSoilStorage])
        if Proc == 'Nash':
            self.routingSystem=LinearChannel(pars=[pars[2],pars[3]])
        if Proc == 'UH':
            self.routingSystem=LinearChannel(pars=pars[2],Proc='UH')
        self.Precipitation=np.array(Boundaries[:,0],dtype='float')
        self.EVP=np.array(Boundaries[:,1],dtype='float')
        self.EVR1=np.array([0]*len(self.Precipitation),dtype='float')
        self.EVR2=np.array([0]*len(self.Precipitation),dtype='float')
        self.SurfaceStorage=np.array([InitialConditions[0]]*(len(self.Precipitation)+1),dtype='float')
        self.SoilStorage=np.array([InitialConditions[1]]*(len(self.Precipitation)+1),dtype='float')
        self.NetRainfall=np.array([0]*len(self.Precipitation),dtype='float')
        self.Infiltration=np.array([0]*len(self.Precipitation),dtype='float')
        self.Runoff=np.array([0]*len(self.Precipitation),dtype='float')
        self.Q=np.array([0]*len(self.Precipitation),dtype='float')
    def computeRunoff(self): 
        j=0
        indexes=list()
        for row in list(self.Precipitation):
            if(self.Precipitation[j]!=0):
                indexes.append(j)
            else:
                if(len(indexes)>0): #Activa Rutina de cómputo SCS
                    print("ponding")
                    self.soilSystem.Precipitation=self.Precipitation[min(indexes):max(indexes)+1]
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
                    indexes=list()
                if(len(indexes)==0): #Activa rutina Cómputo EVR y realiza balance en reservorio de abstracción superficial y reservorio de retención de agua en el suelo
                    print("drying")
                    self.EVR1[j]=min(self.SurfaceStorage[j]/self.maxSurfaceStorage*self.EVP[j],self.SurfaceStorage[j])
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

class HOSH4P2L:
    """
    Modelo Operacional de Transformación de Precipitación en Escorrentía de 4/6 parámetros (estimables), con 2 capas de suelo. Hidrología Operativa Síntesis de Hidrograma. Método NRCS, perfil de suelo con 2 reservorios de retención (zona superior) y un reservorio linear (zona inferior). Rutea utilizando una función respuesta de pulso unitario arbitraria o mediante na cascada de Nash (se debe especificar tiempo de residencia y número de reservorios)
    """
    type='PQ Model'
    def __init__(self,pars,Boundaries=[0],InitialConditions=[[0],[0]],Proc='Nash'):
        self.RoutingProc=Proc
        self.maxSurfaceStorage=pars[0]
        self.maxSoilStorage=pars[1]
        self.soilSystem=SCSReservoirs(pars=[self.maxSurfaceStorage,self.maxSoilStorage])
        self.phi=pars[2]
        self.kb=pars[3]
        if self.RoutingProc == 'Nash':
            self.tr=pars[4]
            self.n=pars[5]
        if self.RoutingProc == 'UH':
            self.u=pars[4]
        self.Precipitation=np.array(Boundaries[:,0],dtype='float')
        self.EVP=np.array(Boundaries[:,1],dtype='float')
        self.EVR1=np.array([0]*len(self.Precipitation),dtype='float')
        self.EVR2=np.array([0]*len(self.Precipitation),dtype='float')
        self.SurfaceStorage=np.array([InitialConditions[0]]*(len(self.Precipitation)+1),dtype='float')
        self.SoilStorage=np.array([InitialConditions[1]]*(len(self.Precipitation)+1),dtype='float')
        self.NetRainfall=np.array([0]*len(self.Precipitation),dtype='float')
        self.Infiltration=np.array([0]*len(self.Precipitation),dtype='float')
        self.Runoff=np.array([0]*len(self.Precipitation),dtype='float')
        self.Q=np.array([0]*len(self.Precipitation),dtype='float')
    def computeRunoff(self): 
        j=0
        indexes=list()
        for row in list(self.Precipitation):
            if(self.Precipitation[j]!=0):
                indexes.append(j)
            else:
                if(len(indexes)>0): #Activa Rutina de cómputo SCS
                    print("ponding")
                    self.soilSystem.Precipitation=self.Precipitation[min(indexes):max(indexes)+1]
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
                    indexes=list()
                if(len(indexes)==0): #Activa rutina Cómputo EVR y realiza balance en reservorio de abstracción superficial y reservorio de retención de agua en el suelo
                    print("drying")
                    self.EVR1[j]=min(self.SurfaceStorage[j]/self.maxSurfaceStorage*self.EVP[j],self.SurfaceStorage[j])
                    self.NetRainfall[j]=max(0,self.Precipitation[j]-self.EVR1[j]+self.SurfaceStorage[j]-self.maxSurfaceStorage)
                    self.EVR2[j]=computeEVR(self.NetRainfall[j],self.EVP[j]-self.EVR1[j],self.SoilStorage[j],self.maxSoilStorage)
                    self.SurfaceStorage[j+1]=waterBalance(self.SurfaceStorage[j],self.Precipitation[j],self.EVR1[j]+self.NetRainfall[j])
                    self.Runoff[j]=max(0,self.NetRainfall[j]-self.EVR2[j]+self.SoilStorage[j-1]-self.maxSoilStorage)
                    self.SoilStorage[j+1]=waterBalance(self.SoilStorage[j],self.NetRainfall[j],self.EVR2[j]+self.Runoff[j])
            j=j+1            
    def computeOutFlow(self):
        if self.RoutingProc == 'Nash':
            self.routingSystem=LinearChannel(pars=[self.tr,self.n],Boundaries=apportion(self.Runoff,self.phi))
        if self.RoutingProc != 'Nash':
            self.routingSystem=LinearChannel(pars=self.u,Boundaries=apportion(self.Runoff,self.phi),Proc='UH')
        self.routingSystem.computeOutFlow()
        self.groundwaterSystem=LinearReservoirCascade(pars=[self.kb,1],Boundaries=apportion(self.Runoff,1-self.phi))
        self.groundwaterSystem.computeOutFlow()
        self.Q=self.routingSystem.Outflow[0:len(self.Runoff)]+self.groundwaterSystem.Outflow[0:len(self.Runoff)]
    def executeRun(self):
        self.computeRunoff()
        self.computeOutFlow()

class GR4J:
    """
    Modelo Operacional de Transformación de Precipitación en Escorrentía de Ingeniería Rural de 4 parámetros (CEMAGREF). A diferencia de la versión original, la convolución se realiza mediante producto de matrices. Parámetros: Máximo almacenamiento en reservorio de producción, tiempo al pico (hidrograma unitario),máximo alamcenamiento en reservorio de propagación, coeficiente de intercambio.
    """
    type='PQ Model'
    def __init__(self,pars,Boundaries=[0],InitialConditions=[[0],[0]],Proc='CEMAGREF SH'):
        self.prodStoreMaxStorage=pars[0]
        self.T=pars[1]
        self.u1=grXDistribution(self.T,distribution='SH1')
        self.u2=grXDistribution(self.T,distribution='SH2')
        self.routStoreMaxStorage=pars[2]
        if not pars[3]:
            self.waterExchange=0
        else:
            self.waterExchange=pars[3]
        self.Precipitation=np.array(Boundaries[:,0],dtype='float')
        self.EVP=np.array(Boundaries[:,1],dtype='float')
        self.Runoff=np.array([0]*len(self.Precipitation),dtype='float')
        self.Q=np.array([0]*len(self.Precipitation),dtype='float')
        self.prodStore=ProductionStoreGR4J(pars=[self.prodStoreMaxStorage],Boundaries=makeBoundaries(self.Precipitation,self.EVP))
    def computeRunoff(self):
        self.prodStore.computeOutFlow()
        self.Runoff=self.prodStore.Runoff
    def computeOutFlow(self):
        self.channel1=LinearChannel(pars=self.u1,Boundaries=apportion(0.9,self.Runoff),Proc='UH')
        self.channel2=LinearChannel(pars=self.u2,Boundaries=apportion(0.1,self.Runoff),Proc='UH')
        self.channel1.computeOutFlow()
        self.routStore=RoutingStoreGR4J(pars=[self.routStoreMaxStorage,self.waterExchange],Boundaries=self.channel1.Outflow)
        self.routStore.computeOutFlow()
        self.channel2.computeOutFlow()
        j=min(len(self.routStore.Runoff),len(self.channel2.Outflow))
        self.Q=self.routStore.Runoff[0:j]+self.channel2.Outflow[0:j]
    def executeRun(self):
        self.computeRunoff()
        self.computeOutFlow()

if __name__ == "__main__":
    import sys


    
