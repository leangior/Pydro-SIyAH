import getSeriesApp as getSeries
import DataDrivenMethods as dd
import matplotlib.pyplot as plt
serie=getSeries.getSerie(155,"1991-01-01","2025-03-30")
anomaly=dd.getCivilAnom(serie)
analogies=dd.getAnalogiesScores(anomaly)
forecast=dd.getAnalogiesValues(analogies,serie)
results=dd.getCentralTrendandForecasts(forecast,serie)

ts=results['centralTrendBiasAdjusted']
x=ts.index.to_numpy()
y1=ts.values.flatten()+3*results['rmse'] #aprox 10% por Chebyshov k=3
y2=ts.values.flatten()-3*results['rmse']
plt.fill_between(x,y1,y2,color='red',alpha=.40)
y1=ts.values.flatten()+4.5*results['rmse'] #aprox 5% por Chebyshov k=4.5
y2=ts.values.flatten()-4.5*results['rmse']
plt.fill_between(x,y1,y2,color='red',alpha=.25)
y1=ts.values.flatten()+7*results['rmse'] #aprox 2% por Chebyshov k=7
y2=ts.values.flatten()-7*results['rmse']
plt.fill_between(x,y1,y2,color='red',alpha=.15)
plt.plot(results['obsSerie'],'o',color='blue')
plt.plot(results['analogies'])
y1=ts.values.flatten()+9*results['rmse'] #aprox 1% por Chebyshov k=7
y2=ts.values.flatten()-9*results['rmse']
plt.fill_between(x,y1,y2,color='red',alpha=.10)
plt.plot(ts,color='blue',linewidth=1.5,linestyle='dashed')
plt.plot(ts,'o',color='red')
plt.ylabel("H[m]")
plt.xlabel("t")
plt.show()

Le puse como parámetro la fecha de pronóstico, cuestión de facilitar el análisis retrospectivo. Aquí perspectiva para el Paraná en Rosario, si hubiera sido emitida el 30/6. No hubiera estado tan mal, che... Igual tendría que realizar un análisis más exhaustivo y me inclino por modelos incluyendo más regresores (ya que hay más hidrómetros y cuencas afleuntes en tramo.. este año por que no escupieron ni una gota de ferné...). Abz grande y buen fin de semana!

Le puse como parámetro la fecha de pronóstico, cuestiónd e facilitar el análisis retrospectivo. Aquí perspectiva para el Uruguay en San javier 
