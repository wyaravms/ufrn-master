
# Valores para os parâmetros
p<-0.40
pho<-0.30
n<-150
size<-100

beta1<-p*(1-pho)
alpha<-beta1+pho

y.aux <- numeric(0)

#Geração de uma sequencia que segue o processo Binomial AR(1)
yn <- function(size,n,alpha,beta1,p){
  y.aux[1] <- rbinom(1,n,p)
  for(i in 2:size)
    y.aux[i] <- rbinom(1,y.aux[(i-1)],alpha) + rbinom(1,n-y.aux[(i-1)],beta1)
  y.aux
}

#Série gerada
y <- yn(size,n,alpha,beta1,p)

plot.ts(y,lwd=1,ylab="",xlab="")
points(y, pch=19, cex=0.5, lwd=1) #Destaca cada ponto da série


mean(y)
var(y) 

#Função de autocorrelação e autocorrelação parcial teorica
m <- 20 #defasagem maxima
par(mfrow=c(1,2))
plot(pho^(1:m),type="h", main = "Função de autocorrelação teorica",ylab = "",xlab="Defasagem")
plot(c(p,rep(0,m-1)),type="h", main = "Função de autocorrelação parcial teorica",ylab = "",xlab="Defasagem")

library("forecast")
#Função de autocorrelação e autocorrelação parcial amostral
par(mfrow=c(1,2))
#tsdisplay(y)

acf(y,lag.max = m)
pacf(y, lag.max = m)

#Yule-Walker
p.YW<-sum(y)/(n*(length(y)+1))
pp <- acf(y,type="partial",plot=FALSE,lag.max=1)
pho.YW <- pp$acf
c(p.YW,pho.YW)

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
c(p.MQC,pho.MQC)


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
c(p.MVC,pho.MVC)

#Parâmetros estimados
c(p,pho)
YW <- round(c(p.YW,pho.YW),3)
MQC <- round(c(p.MQC,pho.MQC),3)
MVC <- round(c(p.MVC,pho.MVC),3)
rbind(c(p,pho),YW,MQC,MVC)
