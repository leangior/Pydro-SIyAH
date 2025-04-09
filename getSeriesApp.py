from a5client import a5 
from a5client.a5 import Crud
import a5_client_utils 
import datetime
import yaml
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
from typing import Optional, Union, List, Tuple

#Load Config File (Services and Credentials) --> Credential must be stated in config.py from a5client
# def loadConfig(configFile="config/config.yml"):
#     with open(configFile) as file:
#         service=yaml.safe_load(file)
#         return(service)

#Get Serie by Id and Type
def getSerie(serieId : int,timeStart,timeEnd,aggStamp : str='D',configFile: str="config/config.yml" ,serieType: str='puntual' ) -> pd.DataFrame:
    # service=loadConfig(configFile)
    client=Crud()
    serie=client.readSerie(serieId,timeStart,timeEnd,serieType)
    colName=serie['estacion']['nombre']+"_"+serie['var']['var']
    serie=a5.observacionesListToDataFrame(serie["observaciones"])
    serie.index=serie.index.tz_convert(None)
    serie=serie.resample(aggStamp).mean()
    serie.columns=[colName]
    return(serie)

def getSeriesDataFrame(seriesId : List[str],timeStart,timeEnd,aggStamp : str='D' ,configFile: str="config/config.yml" ,seriesTypes: List[str]=['puntual']) -> pd.DataFrame:
    # service=loadConfig(configFile)
    client=Crud()
    s=[]
    for i in range(0,len(seriesId)):
        v=client.readSerie(seriesId[i],timeStart,timeEnd,seriesTypes[i])
        colName=v['estacion']['nombre']+"_"+v['var']['var']
        v=a5.observacionesListToDataFrame(v["observaciones"])
        v.index=v.index.tz_convert(None)
        v=v.resample(aggStamp).mean()
        v.columns=[colName]
        s.append(v)
    s=pd.concat(s,axis=1)
    return(s)