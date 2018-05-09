
rm(list=ls(all=TRUE))
set.seed(7)
setwd("C:\\Users\\WYARAVMS\\Google Drive")

# Função da distribuição posterior para os parâmetros
  posterior.burr<-function(dados,int,theta){
       n=length(dados)
       burnin=200
       thin=10
       lpost=function(dados,p,b){
            somalogdados=sum(log(dados))
            somalogdados1=sum(log(1+(dados^b)))
            logpost=n*log(p*b)+(b-1)*(somalogdados)-(p+1)*(somalogdados1)
            logpost=logpost+(0.01-1)*log(p)-(p*0.01)+(0.01-1)*log(b)-(b*0.01)
            logpost
       }
       pmc=array(0,c(burnin+int,1))
       bmc=array(0,c(burnin+int,1))
       pmc[1]= theta[1]
       bmc[1]= theta[2]
       Vp=0.001
       Vb=0.001
       for (i in 2:burnin){
           pest=rgamma(1,pmc[i-1]^2/Vp,pmc[i-1]/Vp)
           best=rgamma(1,bmc[i-1]^2/Vb,bmc[i-1]/Vb)
           alpha=exp(lpost(dados,pest,best)-lpost(dados,pmc[i-1],bmc[i-1]))
           alpha=alpha/(dgamma(pest,pmc[i-1]^2/Vp,pmc[i-1]/Vp))
           alpha=alpha/(dgamma(best,bmc[i-1]^2/Vb,bmc[i-1]/Vb))
           alpha=alpha*dgamma(pmc[i-1],pest^2/Vp,pest/Vp)
           alpha=alpha*dgamma(bmc[i-1],best^2/Vb,best/Vb)
           u=runif(1)
           if (u<alpha){
               pmc[i]=pest
               bmc[i]=best
           }
            else{
               pmc[i]=pmc[i-1]
               bmc[i]=bmc[i-1]
               }
               #if ((i%%100)==0)
               #print(i/(burnin+thin*int))
       }
       pmcb=array(0,c(int))
       bmcb=array(0,c(int))
       j=1
          for (i in (burnin+1):(burnin+thin*int)){
              pest=rgamma(1,pmc[i-1]^2/Vp,pmc[i-1]/Vp)
              best=rgamma(1,bmc[i-1]^2/Vb,bmc[i-1]/Vb)
              alpha=exp(lpost(dados,pest,best)-lpost(dados,pmc[i-1],bmc[i-1]))
              alpha=alpha/(dgamma(pest,pmc[i-1]^2/Vp,pmc[i-1]/Vp))
              alpha=alpha/(dgamma(best,bmc[i-1]^2/Vb,bmc[i-1]/Vb))
              alpha=alpha*dgamma(pmc[i-1],pest^2/Vp,pest/Vp)
              alpha=alpha*dgamma(bmc[i-1],best^2/Vb,best/Vb)
              #if (is.nan(alpha)) {
              #alpha=0
              #}
              u=runif(1)
              if (u<alpha){
                 pmc[i]=pest
                 bmc[i]=best
              }
                 else{
                      pmc[i]=pmc[i-1]
                      bmc[i]=bmc[i-1]
                 }
              if ((i%%thin)==0){
                 pmcb[j]=pmc[i]
                 bmcb[j]=bmc[i]
                 j=j+1
              }
              #if ((i%%100)==0)
        #print(i/(burnin+thin*int))
		}
        cbind(pmcb,bmcb)
  }

# Tamanho da Amostra
n<-25
# Parâmetros
p<-5
b<-5
#Número de Réplicas
R<-3000
#Número de iterações Boostrap
B<-500


# Resultados para os estimadores Bayesianos
result1=matrix(rep(NA,2*R),ncol=2)
result1b=matrix(rep(NA,2*B),ncol=2)
thetaest1=matrix(rep(NA,2*R),ncol=2)

# Resultados para os estimadores de Maxima Verossimilhança
result2=matrix(rep(NA,2*R),ncol=2)
result2b=matrix(rep(NA,2*B),ncol=2)
thetaest2=matrix(rep(NA,2*R),ncol=2)

norm.lik<-function(theta,x){
  p<-theta[1]
  b<-theta[2]
  n<-length(x)
  logl<-n*log(p*b)+(b-1)*sum(log(x))  -(p+1)*sum(log(1+x^b))
  return(-logl)
  }

i=1
  while(i <= R){
    u<-runif(n)
    p<-5
    b<-5
    x<-(((1-u)^(-1/p))-1)^(1/b)

  Max=try(optim(par=c(p,b),fn=norm.lik,x=x,method="BFGS"),T)

  if(class(Max)!="try-error"){
    thetaMax=Max$par
    result2[i,]=thetaMax

    post=posterior.burr(x,1000, thetaMax)
    thetaPost=c(mean(post[,1]), mean(post[,2]))
    result1[i,]=thetaPost

# Início do Boostrap para cada Réplica de Monte Carlo
j=1
   while(j<=B){
      ub=runif(n)
      pb=thetaMax[1]
      bb=thetaMax[2]
      xb<-(((1-ub)^(-1/pb))-1)^(1/bb)

      Maxb=try(optim(par=c(pb,bb),fn=norm.lik,x=xb,method="BFGS"),T)
      if(class(Maxb)!="try-error"){
        thetaMaxb=Maxb$par
        result2b[j,]=thetaMaxb

        postb=posterior.burr(x,1000,thetaMaxb)
        thetaPostb=c(mean(postb[,1]), mean(postb[,2]))
        result1b[j,]=thetaPostb

        j=j+1
        print(j)
      }
   }

	Media1=c(mean(result1b[,1]),mean(result1b[,2])) #bayes
	Media2=c(mean(result2b[,1]),mean(result2b[,2])) #MV

	bias1=Media1-thetaPost        #bayes
	thetaest1[i,]=thetaPost-bias1

	bias2=Media2-thetaMax         #MV
	thetaest2[i,]=thetaMax-bias2
    i=i+1
  print(j)
  }
  print(i)
}

library(xtable)

sink('resultados para os estimadores corrigidos e réplicas - n=25-p=1-1.txt')

# Resultados do Viés Corrigido por Boostrap
Mediacorr1=c(mean(thetaest1[,1]),mean(thetaest1[,2]))
Varcorr1=c(var(thetaest1[,1]),var(thetaest1[,2]))
theta1=c(p,b)
biascorr1=Mediacorr1-theta1
relbiascorr1=biascorr1/theta1*100
EQMcorr1=biascorr1^2+Varcorr1

Mediacorr2=c(mean(thetaest2[,1]),mean(thetaest2[,2]))
Varcorr2=c(var(thetaest2[,1]),var(thetaest2[,2]))
theta2=c(p,b)
biascorr2=Mediacorr2-theta2
relbiascorr2=biascorr2/theta2*100
EQMcorr2=biascorr2^2+Varcorr2


print('Bayesian Estimatot')
mresultscorr1=matrix(c(Mediacorr1,biascorr1,Varcorr1,EQMcorr1),2,4)
dimnames(mresultscorr1)=list(NULL,c("Mediacorr1","biascorr1","Variancia1","EQMcorr1"))
xtable(mresultscorr1,digits = 4)

print('Maximum Likelihood')
mresultscorr2=matrix(c(Mediacorr2,biascorr2,Varcorr2,EQMcorr2),2,4)
dimnames(mresultscorr2)=list(NULL,c("Mediacorr2","biascorr2","Varcorr2","EQMcorr2"))
xtable(mresultscorr2,digits = 4)


# Resultados das Réplicas de Monte Carlo sem correção
Media1=c(mean(result1[,1]),mean(result1[,2]))
Var1=c(var(result1[,1]),var(result1[,2]))
theta1=c(p,b)
bias1=Media1-theta1
relbias1=bias1/theta1*100
EQM1=bias1^2+Var1

Media2=c(mean(result2[,1]),mean(result2[,2]))
Var2=c(var(result2[,1]),var(result2[,2]))
theta2=c(p,b)
bias2=Media2-theta2
relbias2=bias2/theta2*100
EQM2=bias2^2+Var2

print('Bayesian Estimator')
mresults1=matrix(c(Media1,bias1,Var1,EQM1),2,4)
dimnames(mresults1)=list(NULL,c("Media1","bias1","Var1","EQM1"))
xtable(mresults1,digits = 4)


print('Maximum Likelihood')
mresults2=matrix(c(Media2,bias2,Var2,EQM2),2,4)
dimnames(mresults2)=list(NULL,c("Media2","bias2","Var2","EQM2"))
xtable(mresults2,digits = 4)

sink()
