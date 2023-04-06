#/usr/bin/python3
#Librería de métodos para evaluación de hindcast, 2023
import math
import numpy  as np
import pandas as pd
import scipy.stats


#Cómputo de Error Porcentual Absoluto Medio. Toma por entrada el valor observado (flotante) y un ensamble de pronóstico (lista o array 1D). Computa el error porcentual absoluto de cada miembro de ensamble en relación al valor observado (ape). Finalmente computa su media aritmética (mape). Es un indicador del error de volumen.     
def mape(obs,ensembles=[]):
    ensembles=np.array(ensembles,dtype='float')
    obs=np.array([obs]*len(ensembles),dtype='float')
    ape=abs((ensembles-obs)/obs)*100
    mape=sum(ape)/len(ape)
    return(mape)

# Old School way
# def mape(obs,ensembles=[]):
#     mape=0
#     for member in ensembles:
#         mape=mape+abs((member-obs)/obs)*100
#     return(mape/len(ensembles))

#Cómputo de Rango Percentil. Toma por entrada el valor observado (flotante) y un ensamble de pronóstico (lista o array 1D). Computa la posición del valor observado en la lista de miembros de ensamble. Computa el valor de frecuencia (acumulada) en la distribución empírica observada en el ensamble. Es un indicador de la 'calidad' del ensamble (cuán hábil es para simular la variabilidad observable), esto es: si el ensamble subestima la variabilidad o si la sobreestima. Puede interpretarse como un indicador de tendencia (sesgo). El paquete scipy ofrece una alternativa similar.
def prs(obs,ensembles=[]):
    prs=scipy.stats.percentileofscore(np.array(ensembles),obs)
    return(prs)

# Old School way
# def prs(obs,ensembles=[]):
#     rank=0
#     for member in ensembles:
#         if(member < obs):
#             rank=rank+1
#     return(rank/len(ensembles))

#Cómputo del coeficiente de asociación no paramétrico 'Tau' de Kendall. Se adopta la función del paquete scipy, para cualquier consulta de documentación. Las entradas están constituídas por una lista de valores observados y otra lista de valores simulados. Es un indicador robusto del grado de asociación lineal (dependencia lineal por proporcionalidad o aditividad), con dominio [-1,1] (los valores límites indican asociación perfecta negativa o asociación perfecta positiva, dependencia lineal absoluta, por otro lado 0 indica indepencia lineal). 
def tauKendall(obs=[],sim=[]):
    tau=scipy.stats.kendalltau(sim,obs)[0]
    return(tau)

#Cómputo del coeficiente de skill Nash-Sutcliffe. Indicador de la asociación 1:1 entre valores simulados y observados. 
def nashSutcliffeScore(obs=[],sim=[]):
    obs=np.array(obs)
    sim=np.array(sim)
    ns=1-sum((obs-sim)**2)/sum((obs-np.mean(obs))**2)
    return(ns)