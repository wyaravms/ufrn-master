
setwd("C:\\Users\\Wyara\\Dropbox\\series-temporais-seminario")

p<-0.40
pho<-0.30
n<-50
size<-50

beta1<-p*(1-pho)
alpha<-beta1+pho
y.aux <- numeric(0)
R<-1000

# Resultados para os estimadores 
resultyw=matrix(rep(NA,2*R),ncol=2)
resultmqc=matrix(rep(NA,2*R),ncol=2)
resultmvc=matrix(rep(NA,2*R),ncol=2)

start.time <- Sys.time()
j=1
  while(j <= R){

	#Geração de uma sequencia que segue o processo Binomial AR(1)
	yn <- function(size,n,alpha,beta1,p){
		y.aux[1] <- rbinom(1,n,p)
		for(i in 2:size)
		y.aux[i] <- rbinom(1,y.aux[(i-1)],alpha) + rbinom(1,n-y.aux[(i-1)],beta1)
	y.aux
}

	#Série gerada
	y <- yn(size,n,alpha,beta1,p)
	
	#Yule-Walker
	p.YW<-sum(y)/(n*(length(y)+1))
	pp <- acf(y,type="partial",plot=FALSE,lag.max=1)
	pho.YW <- pp$acf
	est1=c(p.YW,pho.YW)
	resultyw[j,]=est1
	
	#Minimos quadrados condicionais
	p.MQC <- sum(y)/(n*(length(y)))

	a1 <- 0
	b1 <- 0
	for (i in 2:size){   
		a<-(y[i]-n*p.MQC)*(y[i-1]-n*p.MQC)
		a1<-a1+a
		b<-((y[i-1]-n*p.MQC)^2)
		b1<-b1+b
	}

	pho.MQC<-a1/b1
	#pho.MQC <- (sum(y[2:n])-n*p.MQC*sum(y[1:(n-1)]))/(n-1)
	
	est2<-c(p.MQC,pho.MQC)
	resultmqc[j,]=est2
	
	#Máxima verossimilhança condicional
	#------------------------------------

	#Função de probabiidade condicional
	fdp <- function(yt,alpha=alpha,beta1=beta1,yt1,n){
		mini <- min(yt,yt1)
		maxi<-max(0,(yt+yt1-n))
		w <- maxi:mini
		t<-choose(yt1,w)*choose(n-yt1,yt-w)*(alpha^w)*((1-alpha)^(yt1-w))*(beta1^(yt-w))*((1-beta1)^(n-yt1+w-yt))
		sum(t)
	}	

	#Função de log-verossimilhança
	library(maxLik)

	aux <- numeric(0)
	mlog <- function(theta,y,size){
	for(t in 2:size){
		aux[(t-1)] <- log(fdp(yt=y[t],alpha=theta[1],beta1=theta[2],yt1=y[(t-1)],n))
	}
  sum(aux)
}
	
	beta.MQC<-p.MQC*(1-pho.MQC)
	alpha.MQC<-pho.MQC+beta.MQC
	theta.start <- c(alpha.MQC,beta.MQC)
	
	out <- maxLik(mlog,grad = NULL, hess = NULL, start=theta.start,y=y,size=size)  
	alpha.MVC<-out$estimate[1]
	beta.MVC<-out$estimate[2]
	pho.MVC <-alpha.MVC-beta.MVC 
	p.MVC <-beta.MVC/(1-pho.MVC)
	
	est3<-c(p.MVC,pho.MVC)
	resultmvc[j,]=est3
	
	j=j+1
  
  }
  end.time <- Sys.time()

	sink('time to finished- 1-n=50.txt')
	start.time
	end.time

	time.taken <- end.time - start.time
	time.taken

	sink()

  sink('resultados dos estimadores para o Monte Carlo- 1 - n=50.txt')
  
  theta=c(p,pho)
  
  MediaYW=c(mean(resultyw[,1]),mean(resultyw[,2]))
  VarYW=c(var(resultyw[,1]),var(resultyw[,2]))
  biasYW=MediaYW-theta
  EQMYW=biasYW^2+VarYW
  
  MediaMQC=c(mean(resultmqc[,1]),mean(resultmqc[,2]))
  VarMQC=c(var(resultmqc[,1]),var(resultmqc[,2]))
  biasMQC=MediaMQC-theta
  EQMMQC=biasMQC^2+VarMQC
  
  MediaMVC=c(mean(resultmvc[,1]),mean(resultmvc[,2]))
  VarMVC=c(var(resultmvc[,1]),var(resultmvc[,2]))
  biasMVC=MediaMVC-theta
  EQMMVC=biasMVC^2+VarMVC
  
  mresultsYW=matrix(c(MediaYW,VarYW,biasYW,EQMYW),2,4)
  dimnames(mresultsYW)=list(NULL,c("MediaYW","varYW","biasYW","EQMYW"))
  mresultsYW   

  mresultsMQC=matrix(c(MediaMQC,VarMQC,biasMQC,EQMMQC),2,4)
  dimnames(mresultsMQC)=list(NULL,c("MediaMQC","varMQC","biasMQC","EQMMQC"))
  mresultsMQC
  
  mresultsMVC=matrix(c(MediaMVC,VarMVC,biasMVC,EQMMVC),2,4)
  dimnames(mresultsMVC)=list(NULL,c("MediaMVC","varMVC","biasMVC","EQMMVC"))
  mresultsMVC
  
  sink()
  
  
  