lambda <- 2
p <- 0.1 #Fazer grafico com p=0.8, 0.9, 0.95
mu <- lambda/(1-p) 
y.aux <- numeric(0)
n <- 150


#Geração de uma sequencia que segue o processo INAR(1)
yn <- function(n){
  y.aux[1] <- rpois(1,mu)
  for(i in 2:n)
    y.aux[i] <- rbinom(1,y.aux[(i-1)],p) + rpois(1,mu)
  y.aux
}

y <- yn(n) #Série gerada

plot.ts(y,lwd=1,ylab="",xlab="")
points(y, pch=19, cex=0.5, lwd=1) #Destaca cada ponto da série


mean(y)
var(y) #Esperamos que estes dois valores sejam próximos


#Função de autocorrelação e autocorrelação parcial teorica
m <- 20 #defasagem maxima
par(mfrow=c(1,2))
plot(p^(1:m),type="h", main = "Função de autocorrelação teorica",ylab = "",xlab="Defasagem")
plot(c(p,rep(0,m-1)),type="h", main = "Função de autocorrelação parcial teorica",ylab = "",xlab="Defasagem")


#Função de autocorrelação e autocorrelação parcial amostral
par(mfrow=c(1,1))
tsdisplay(y)

acf(y,lag.max = m)
pacf(y, lag.max = m)

#Estimação dos parâmetros
#--------------------------

#Yule-Walker
ybarra <- mean(y)
pp <- acf(y,type="partial",plot=FALSE,lag.max=1)
p.YW <- pp$acf
lamb.YW <- (1-p.YW)*ybarra
c(p.YW,lamb.YW)


#Minimos quadrados condicionais
p.MQC <- (sum(y[2:n]*y[1:(n-1)])-sum(y[2:n])*sum(y[1:(n-1)])/(n-1))/(sum(y[1:(n-1)]^2)-(sum(y[1:(n-1)]))^2/(n-1))
lamb.MQC <- (sum(y[2:n])-p.MQC*sum(y[1:(n-1)]))/(n-1)
c(p.MQC,lamb.MQC)

#Máxima verossimilhança condicional
#------------------------------------

#Função de probabiidade condicional
fdp <- function(yt,alpha=alph,lamb=lambd,yt1){
  mini <- min(yt,yt1)
  w <- 0:mini
  t <- lamb^(yt-w)*choose(yt1,w)*alpha^w*(1-alpha)^(yt1-w)/(factorial(yt-w))
  exp(-lamb)*sum(t)
}

#Função de log-verossimilhança
library(maxLik)

aux <- numeric(0)
mlog <- function(theta,y,n){
  for(t in 2:n){
    aux[(t-1)] <- log(fdp(yt=y[t],alpha=theta[1],lamb=theta[2],yt1=y[(t-1)]))
  }
  sum(aux)
}

theta.start <- c(p.MQC,lamb.MQC)
out <- maxLik(mlog,grad = NULL, hess = NULL, start=theta.start,y=y,n=n)
p.MVC <- out$estimate[1]
lamb.MVC <- out$estimate[2]
c(p.MVC,lamb.MVC)

#Parâmetros estimados
c(p,lambda)
YW <- round(c(p.YW,lamb.YW),3)
MQC <- round(c(p.MQC,lamb.MQC),3)
MVC <- round(c(p.MVC,lamb.MVC),3)
rbind(c(p,lambda),YW,MQC,MVC)

#-------------------------------------------
#----------------APLICAÇÂO------------------
#-------------------------------------------


#--------------------------------------------------------------------------------------
#                                  DADOS POLIOMELITIS
#--------------------------------------------------------------------------------------
polio <-c (0,  1,  0,  0,  1,	3,	9,	2,	3,	5,	3,	5, 
           2,	2,	0,	1,	0,	1,	3,	3,	2,	1,	1,	5, 
           0,	3,	1,	0,	1,	4,	0,	0,	1,	6,	14,	1, 
           1,	0,	0,	1,	1,	1,	1,	0,	1,	0,	1,	0, 
           1,	0,	1,	0,	1,	0,	1,	0,	1,	0,	0,	2, 
           0,	1,	0,	1,	0,	0,	1,	2,	0,	0,	1,	2, 
           0,	3,	1,	1,	0,	2,	0,	4,	0,	2,	1,	1, 
           1,	1,	0,	1,	1,	0,	2,	1,	3,	1,	2,	4, 
           0,	0,	0,	1,	0,	1,	0,	2,	2,	4,	2,	3, 
           3,	0,	0,	2,	7,	8,	2,	4,	1,	1,	2,	4, 
           0,	1,	1,	1,	3,	0,	0,	0,	0,	1,	0,	1, 
           1,	0,	0,	0,	0,	0,	1,	2,	0,	2,	0,	0, 
           0,	1,	0,	1,	0,	1,	0,	2,	0,	0,	1,	2, 
           0,	1,	0,	0,	0,	1,	2,	1,	0,	1,	3,	6) 


strike = c(	1,	2,	3,	5,	6 ,	11,	9,	11,	14,	9 ,	6,	4,
            4,	4,	7,	7,	4 , 4 , 5,	8 , 9 , 10,	5,	3,
            4,	9,	7,	7,	4 , 11,	7,	5 , 6 , 8 , 7,	5,
            6,	4,	4,	7,	12,	3 , 5,	3 , 6 , 5 , 3,	2,
            1,	3,	2,	0,	4 , 6 , 6,	7 , 4 , 7 , 7,	6,
            6,	5,	3,	4,	6 , 6 , 6,	3 , 5 , 2 , 2,	1,
            1,	2,	4,	7,	4 , 8 , 6,	8 , 10,	12,	3,	3,
            2,	1,	4,	5,	8 , 5 , 3,	4 , 3 , 4 , 1,	2,
            1,	2,	1,	3,	5 , 3 , 4,	3 , 3 , 3 , 2,	1)

soft <- c( 9,  6,   6,   7,  10,   8,  14,   8,   7,  10,  10,  12,   8,   8,   8,
           8, 13, 12,  14, 13,  13,   8,  13,  10,  12,  12,   9,   8,  13,   9,   8, 
           6,  7,  10,  17,  11,  13,  10,   9,  15,  13,  12,   8,   8,   9,   9,  12,
           9,   5,   9,  10,   6,   8,  17,  16,  17,  16,   8,  10,   7,   8,   7, 
           4,  5,   4,   4,  10,   9,  12,  12,  11,   9,   8,   9,   8,   6,   8,  13,
           13,  10,   7,  17,  14,  10,  12,   6,   4,   7,   8,  10,  16,  15,  10, 
           14, 16,  12,  10,  11,  10,   8,   9,  10,  13,   6,   8,   9,   6,   9,  12,
           8,   9,   5,   6,   9,   9,  13,  12,  10,   9,   7)


dislocation <- c(0, 0, 1, 1, 1, 1, 1, 1, 0, 1, 1, 2, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0,
                 0, 0, 0, 0, 2, 1, 0, 0, 2, 2, 2, 2, 1, 0, 0, 1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 1, 2, 2, 1,
                 1, 0, 1, 3, 1, 1, 1, 1, 2, 1, 1, 2, 2, 1, 1, 3, 4, 3, 2, 1, 1, 0, 0, 0, 1, 2, 0, 0, 0, 0,
                 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 1, 0, 0, 1, 2, 3, 2)


Weib2007<- c (0, 0, 2, 0, 3, 1, 0, 2, 0, 0, 2, 2, 3, 0, 0, 0, 0, 1, 0, 2, 3, 1, 0, 0, 1, 1, 3, 3, 2, 1, 1, 1, 1, 2, 1, 1, 2,
              0, 1, 1, 2, 2, 2, 2, 2, 1, 1, 1, 3, 3, 1, 0, 0, 2, 3, 1, 0, 1, 0, 0, 1, 0, 0, 1, 1, 1, 1, 1, 0, 2, 2, 4, 2, 2,
              1, 1, 1, 1, 0, 1, 4, 2, 3, 1, 1, 3, 4, 2, 3, 0, 0, 0, 1, 0, 0, 1, 2, 1, 1, 2, 0, 1, 0, 0, 0, 0, 1, 2, 1, 0, 3,
              1, 2, 0, 1, 0, 2, 1, 1, 0, 0, 0, 2, 0, 3, 0, 1, 1, 2, 0, 0, 2, 2, 0, 1, 2, 3, 3, 2, 3, 2, 2, 2, 2, 1, 1, 1, 0,
              2, 3, 2, 4, 1, 2, 4, 2, 2, 2, 0, 2, 2, 1, 1, 1, 0, 1, 0, 1, 0, 0, 1, 1, 1, 2, 3, 5, 2, 2, 1, 4, 5, 3, 0, 0, 2,
              1, 2, 0, 2, 2, 2, 1, 2, 2, 1, 1, 0, 1, 2, 2, 1, 0, 3, 2, 2, 2, 2, 3, 0, 1, 1, 0, 1, 3, 3, 1, 1, 2, 1, 0, 1, 0,
              1, 8, 0, 0, 1, 3, 0, 0, 2, 1, 1, 0, 1, 1, 0, 1, 1, 1, 0)  	  	 	   

n1 <- length(polio)
n2 <- length(strike)
n3 <- length(soft)
n4 <- length(dislocation)
n5 <- length(Weib2007)

POLIO <- c(min(polio),max(polio),median(polio),mean(polio),(n1-1)*var(polio)/n1)
STRIKE <- c(min(strike),max(strike),median(strike),mean(strike),(n2-1)*var(strike)/n2)
SOFT <- c(min(soft),max(soft),median(soft),mean(soft),(n3-1)*var(soft)/n3)
DISLOCATION <- c(min(dislocation),max(dislocation),median(dislocation),mean(dislocation),(n4-1)*var(dislocation)/n4)
WEIB <- c(min(Weib2007),max(Weib2007),median(Weib2007),mean(Weib2007),(n5-1)*var(Weib2007)/n5)

medidas <- rbind(POLIO,STRIKE,SOFT,DISLOCATION,WEIB)
colnames(medidas) <- c("min","max","mediana","média","var")
print(medidas,digits=3, include.rownames=T)

library(forecast)
y2 <- dislocation

tsdisplay(y2,main="STRIKE")

theta.start=c(0.1,1)
out2 <- maxLik(mlog,grad = NULL, hess = NULL, start=theta.start,y=y2,n=n4)
p.MVC.d <- out2$estimate[1]
lamb.MVC.d <- out2$estimate[2]
c(p.MVC.d,lamb.MVC.d)

#EVMC com chute inicial os estimadores MQC

p.MQC <- (sum(y2[2:n4]*y2[1:(n4-1)])-sum(y2[2:n4])*sum(y2[1:(n4-1)])/(n4-1))/(sum(y2[1:(n4-1)]^2)-(sum(y2[1:(n4-1)]))^2/(n4-1))
lamb.MQC <- (sum(y2[2:n4])-p.MQC*sum(y2[1:(n4-1)]))/(n4-1)

theta.start <- c(p.MQC,lamb.MQC)
out2 <- maxLik(mlog,grad = NULL, hess = NULL, start=theta.start,y=y2,n=n4)
p.MVC.d <- out2$estimate[1]
lamb.MVC.d <- out2$estimate[2]
c(p.MVC.d,lamb.MVC.d)

#Residuos

residuos.dis <- numeric(0)
for(i in 2:n4){
  residuos.dis[(i-1)] <- y2[i]-p.MVC.d*y2[(i-1)]-lamb.MVC.d 
}

acf(residuos.dis) #Autocorrelação dos residuos 

mean(residuos.dis)
var(residuos.dis) #lamb.MVC.d*(1+p.MVC.d)









