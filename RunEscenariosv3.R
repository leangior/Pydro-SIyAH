source("/home/maxi/01-PROCEDIMIENTOS/Rsuite/LinearSystems.R")
#Modelación
InitForecast=Sys.Date() #as.Date("2023-05-30")
End=Sys.Date()+30 #as.Date("2023-06-29")
TrazaYaci="TrazaYaciOps"
TrazaFor="TrazaForOps"
type='max'
pos=1 #max=1
#carga trazas
yacireta=read.table(TrazaYaci,head=F,colClasses = c('Date','numeric','numeric','numeric'),fill=T)
colnames(yacireta)=c('t','MAX','MED','MIN')
yaci=xts(na.approx(yacireta$MIN),order.by=as.Date(yacireta$t))
yaci=apply.daily(yaci,FUN=mean)
formosa=read.table(TrazaFor,head=F,colClasses = c('Date','numeric','numeric','numeric'),fill=T)
colnames(formosa)=c('t','MAX','MED','MIN')
#Ejecuta con traza mínima
formo=xts(na.approx(formosa$MIN),order.by=as.Date(formosa$t))
formo=apply.daily(formo,FUN=mean)
min=TransitoParana("2020-01-01",InitForecast,Ahead=length(yaci)-1,Borde_Yacireta = yaci,Borde_Formosa = formo)
min=GetFitsParana(min,"2020-04-01","2022-12-30",End)
#Ejecuta con traza media
yaci=xts(na.approx(yacireta$MED),order.by=as.Date(yacireta$t))
yaci=apply.daily(yaci,FUN=mean)
formo=xts(na.approx(formosa$MED),order.by=as.Date(formosa$t))
formo=apply.daily(formo,FUN=mean)
med=TransitoParana("2020-01-01",InitForecast,Ahead=length(yaci)-1,Borde_Yacireta = yaci,Borde_Formosa = formo)
med=GetFitsParana(med,"2020-04-01","2022-12-30",End)
#Ejecuta con traza max
yaci=xts(na.approx(yacireta$MAX),order.by=as.Date(yacireta$t))
yaci=apply.daily(yaci,FUN=mean)
formo=xts(na.approx(formosa$MAX),order.by=as.Date(formosa$t))
formo=apply.daily(formo,FUN=mean)
max=TransitoParana("2020-01-01",InitForecast,Ahead=length(yaci)-1,Borde_Yacireta = yaci,Borde_Formosa = formo)
max=GetFitsParana(max,"2020-04-01","2022-12-30",End)
#Configuración y exportación de archivos
Corrientes=data.frame(t=index(med$Corrientes),min=min$Corrientes$min,med=med$Corrientes$sim,max=max$Corrientes$max,obs=med$Corrientes$obs)
write.csv(Corrientes,file = "Corrientes_Esc.csv",quote=F,row.names = F)
Barranqueras=data.frame(t=index(med$Barranqueras),min=min$Barranqueras$min,med=med$Barranqueras$sim,max=max$Barranqueras$max,obs=med$Barranqueras$obs)
write.csv(Barranqueras,file = "Barranqueras_Esc.csv",quote=F,row.names = F)
Goya=data.frame(t=index(med$Goya),min=min$Goya$min,med=med$Goya$sim,max=max$Goya$max,obs=med$Goya$obs)
write.csv(Goya,file = "Goya_Esc.csv",quote=F,row.names = F)
Reconquista=data.frame(t=index(med$Reconquista),min=min$Reconquista$min,med=med$Reconquista$sim,max=max$Reconquista$max,obs=med$Reconquista$obs)
write.csv(Reconquista,file = "Reconquista_Esc.csv",quote=F,row.names = F)
LaPaz=data.frame(t=index(med$LaPaz),min=min$LaPaz$min,med=med$LaPaz$sim,max=max$LaPaz$max,obs=med$LaPaz$obs)
write.csv(LaPaz,file = "LaPaz_Esc.csv",quote=F,row.names = F)
Parana=data.frame(t=index(med$Parana),min=min$Parana$min,med=med$Parana$sim,max=max$Parana$max,obs=med$Parana$obs)
write.csv(Parana,file = "Parana_Esc.csv",quote=F,row.names = F)
SantaFe=data.frame(t=index(med$SantaFe),min=min$SantaFe$min,med=med$SantaFe$sim,max=max$SantaFe$max,obs=med$SantaFe$obs)
write.csv(Parana,file = "SantaFe_Esc.csv",quote=F,row.names = F)
Rosario=data.frame(t=index(med$Rosario),min=min$Rosario$min,med=med$Rosario$sim,max=max$Rosario$max,obs=med$Rosario$obs)
write.csv(Parana,file = "SantaFe_Esc.csv",quote=F,row.names = F)
system("zip Escenarios.zip *.csv")
#
LinealFit<-function(serie,Start=Sys.Date()-90,End=Sys.Date()+15){
  x=serie
  Start=as.Date(Start)
  End=as.Date(End)
  Interval=paste0(Start,"/",End)
  sample=xts(serie,order.by=serie$t)[Interval]
  sim=as.numeric(sample$sim)
  obs=as.numeric(sample$obs)
  fit=summary(lm(obs~sim))
  x$fit=fit$coef[1]+fit$coef[2]*x$sim
  x$r.squared=fit$r.squared
  return(x)
}
GetForecast<-function(DataFrame,Start=Sys.Date()-2,End=Sys.Date()+15,Warm=Sys.Date()-90,Stop=End){
  ForecastFrame=LinealFit(DataFrame,Warm,Stop)
  Start=as.Date(Start)
  End=as.Date(End)
  Interval=paste0(Start,"/",End)
  return(xts(ForecastFrame,order.by=ForecastFrame$t)[Interval])
}
PlotEscenario<-function(data,fit=F){
  if(fit==T){
    data=LinealFit(data)
  }
  plot=ggplot(data,aes(x=t))+geom_line(aes(y=max),col='red')+geom_line(aes(y=min),col='darkorange')+geom_line(aes(y=sim),col='black',lty=2)+geom_line(aes(y=obs),col='blue',lty=2)
  if(fit==T){
    plot+geom_line(aes(y=fit),col='darkviolet',lty=5)
  }
  else{
    plot
  }
}
GetDataForecastLine<-function(data,date=Sys.Date()){
  sample=c()
  date=as.Date(date)
  sample[1]=as.numeric(GetForecast(data)[date]$fit)
  sample[2]=as.numeric(GetForecast(data)[date]$sim)
  sample[3]=as.numeric(GetForecast(data)[date]$min)
  sample[4]=as.numeric(GetForecast(data)[date]$max)
  return(sample)
}
GetLinearWeight<-function(data,date=Sys.Date()){
  date=as.Date(date)
  sims=GetDataForecastLine(data,date)
  obs=as.numeric(GetForecast(data)[date]$obs)
  w=c()
  i=1
  for(value in sims){w[i]=(obs-sims[i])^(-2);i=i+1}
  w=w/sum(w)
  if(sum(w)>=.999||sum(w)<=1){
    return(w)
  }
  else{
    print("error en cómputo sum(w)!=1")
  }
}
ComputeCentralTrend<-function(data,init_date=Sys.Date(),lag=7){
  w=GetLinearWeight(data,init_date)
  PreForecast=GetDataForecastLine(data,init_date+lag)
  forecast=c()
  forecast$date=init_date+lag
  forecast$value=PreForecast%*%w
  forecast$range=c(min(PreForecast),max(PreForecast))
  forecast$init=as.numeric(GetForecast(data)[init_date]$obs)
  return(forecast)
}
getBiasModel<-function(analisis,obs,sim,T,type='max'){
  if(type=='max'){
    sub=subset(analisis,obs>T)  
  }
  else{
    sub=subset(analisis,obs<=T)
  }
  model=summary(lm(sub$obs~sub$sim))
  result=list()
  result$model=model
  result$sub=sub
  return(result)
}
calibTravelTimeReach<-function(upstreamSerieId,downstreamSerieId,lag=15,start=Sys.Date()-365*2,end=Sys.Date(),kIni=0.5,kStop=8,n=2,dk=0.1,type='rising',Agg='Daily'){
  upstream=GetDBSerie(upstreamSerieId,start,end,Agg=Agg)
  downstream=GetDBSerie(downstreamSerieId,start,end,Agg=Agg)
  interval=paste0(Sys.Date()-lag,"/",Sys.Date())
  values=list()
  i=0
  for(k in seq(kIni,kStop,dk)){
    i=i+1
    route=LinearTransit(upstream,k=k,n=n)
    data=cbind(downstream,route)
    if(type=='rising'){
      values$data=subset(data,data$downstream>=min(data$downstream[interval])) 
    }
    else{
      values$data=subset(data,data$downstream<=max(data$downstream[interval]))
    }
    values$k[i]=k
    values$r2[i]=summary(lm(as.numeric(values$data$downstream)~as.numeric(values$data$route)))$r.squared
  }
  values$kOpt=values$k[which.max(values$r2)]
  route=LinearTransit(upstream,k=values$kOpt,n=n)
  values$data=cbind(downstream,route)
  if(type=='rising'){
    values$subset=subset(values$data,data$downstream>=min(data$downstream[interval])) 
  }
  else{
    values$subset=subset(values$data,data$downstream<=min(data$downstream[interval]))
  }
  values$fit=summary(lm(as.numeric(values$data$downstream)~as.numeric(values$data$route)))
  values$Qref=min(data$downstream[interval])
  return(values)
}
EvalPredictedValue<-function(model,serie,lag){
  results=list()
  results$time=Sys.Date()+lag
  results$r2=round(model$r.squared,2)
  results$range=ComputeCentralTrend(serie,lag=lag)$range
  results$value=ComputeCentralTrend(serie,lag=lag)$value
  results$fit=model$coef[1]+model$coef[2]*ComputeCentralTrend(serie,lag=lag)$value
  results$ave=mean(c(results$value,results$fit))
  results$base=ComputeCentralTrend(serie,lag=lag)$init
  return(results)
}
getSerieSimAndObsFromTransits<-function(t,vals_nodo){
  serie=xts(order.by=t,vals_nodo)
  colnames(serie)=c('obs','sim')
  return(serie)
}
getErrorModelAndAdjust<-function(t,obs,sim,st=Sys.Date()-30,ed=Sys.Date(),k=1.68){
  st=as.Date(st)
  ed=as.Date(ed)
  s=getSerieSimAndObsFromTransits(t,cbind(obs,sim))
  colnames(s)=c('obs','sim')
  t=seq(st,ed,by='1 day')
  sub=s[t]
  v0=sub$obs[ed]
  x=xts(order.by=ed,v0)
  model=summary(lm(sub$obs~sub$sim))
  model_error=forecast(auto.arima(sub$obs-sub$sim))
  t=seq(ed+1,ed+length(model_error$mean),by='1 day')
  sub=s[t]
  lo=xts(order.by=t,model_error$lower[,1])
  up=xts(order.by=t,model_error$upper[,1])
  u=xts(order.by=t,model_error$mean)
  errorModeledSerie=xts(order.by=t,cbind(sub$sim,sub$sim+model$sigma*k,sub$sim-model$sigma*k,sub$sim+lo,sub$sim+u,sub$sim+up))
  for(i in seq(2,dim(errorModeledSerie)[2])){
    x=cbind(x,xts(order.by=ed,v0))
  }
  errorModeledSerie=rbind(x,errorModeledSerie)
  colnames(errorModeledSerie)=c('sim_mean','sim_up','sim_do','sim_err_do','sim_err_mean','sim_err_up')
  return(errorModeledSerie)
}
lagAndRoute<-function(upSerie,lag=0,k=0.01,n=1,warmUp=30){
    route=upSerie
    index(route)=index(route)+lag
    route=LinearTransit(route,k,n)
    return(route[warmUp:dim(route)[1]])
}
lagAndRouteBiasAdj<-function(upSerie,downSerie,lag=0,k=0.01,n=1,warmUp=30){
  route=lagAndRoute(upSerie,lag,k,n,warmUp)
  series=cbind(downSerie,route)
  series=series[warmUp:dim(series)[1]]
  colnames(series)=c('obs','sim')
  model=summary(lm(series$obs~series$sim))
  message(paste0("Utilizando ajuste de sesgo por regresión lineal R²=",model$r.squared))
  series$sim=model$coef[1]+model$coef[2]*series$sim
  return(series)
}
Bias=list()
Bias$Cor=getBiasModel(Corrientes,Corrientes$obs,Corrientes$sim,ComputeCentralTrend(Corrientes,lag=7)$range[pos],type)
Bias$Bar=getBiasModel(Barranqueras,Barranqueras$obs,Barranqueras$sim,ComputeCentralTrend(Barranqueras,lag=7)$range[pos],type)
Bias$Goya=getBiasModel(Goya,Goya$obs,Goya$sim,ComputeCentralTrend(Goya,lag=7)$range[pos],type)
Bias$Rec=getBiasModel(Reconquista,Reconquista$obs,Reconquista$sim,ComputeCentralTrend(Reconquista,lag=7)$range[pos],type)
Bias$Paz=getBiasModel(LaPaz,LaPaz$obs,LaPaz$sim,ComputeCentralTrend(LaPaz,lag=7)$range[pos],type)
Bias$Par=getBiasModel(Parana,Parana$obs,Parana$sim,ComputeCentralTrend(Parana,lag=7)$range[pos],type)
Bias$Sfe=getBiasModel(SantaFe,SantaFe$obs,SantaFe$sim,ComputeCentralTrend(SantaFe,lag=7)$range[pos],type)
Bias$Ros=getBiasModel(Rosario,Rosario$obs,Rosario$sim,ComputeCentralTrend(Rosario,lag=7)$range[pos],type)
