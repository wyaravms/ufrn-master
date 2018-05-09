
library(VGAM)

posterior.simplex<-function(dados,int){
     thin=10
     burnin=int*thin/2 
     n=length(dados)
     lpost=function(dados,mup,lambdap){
           sumdados=sum(((dados-mup)^2)/(dados*(1-dados)*(mup^2)*(1-mup)^2))
           logpost=(n/2)*log(1/lambdap^2)-((1/(2*lambdap^2))*sumdados)
           logpost=logpost+(0.01-1)*log(mup)+(0.01-1)*log(1-mup)+(0.001-1)*log(lambdap^2)-(0.001*lambdap^2)
           logpost
     }
     mumc=array(0,c(burnin+int,1))
     lamdamc=array(0,c(burnin+int,1))
  
     #usando o chute inicial pelo emv do lambda
     mumc[1]=mean(dados)
     lamdamc[1]=2
  
     Vmu=0.001  
     Vlamb=0.01

     for (i in 2:burnin){
         muest=rbeta(1,(((mumc[i-1]^2)*(1-mumc[i-1]))/Vmu)+((mumc[i-1])/(1-mumc[i-1]))*(mumc[i-1]-1),
         ((((mumc[i-1])*((mumc[i-1])^2))/Vmu)-1+mumc[i-1]))   
         lambest=rgamma(1,lamdamc[i-1]^2/Vlamb,lamdamc[i-1]/Vlamb)
    
         alpha=exp(lpost(dados,muest,lambest)-lpost(dados,mumc[i-1],lamdamc[i-1]))
         alpha=alpha/(dbeta(muest,(((mumc[i-1]^2)*(1-mumc[i-1]))/Vmu)+((mumc[i-1])/(1-mumc[i-1]))*(mumc[i-1]-1),
         ((((mumc[i-1])*((mumc[i-1])^2))/Vmu)-1+mumc[i-1])))
         alpha=alpha/(dgamma(lambest,lamdamc[i-1]^2/Vlamb,lamdamc[i-1]/Vlamb))
         alpha=alpha*dbeta(mumc[i-1],((((muest^2)*(1-muest))/Vmu)+(muest/(1-muest))*(muest-1)),
         (((muest*(muest^2))/Vmu)-1+muest))
         alpha=alpha*dgamma(lamdamc[i-1],lambest^2/Vlamb,lambest/Vlamb)

         u=runif(1)
         if (u<alpha){
            mumc[i]=muest
            lamdamc[i]=lambest
         } else{
            mumc[i]=mumc[i-1]
            lamdamc[i]=lamdamc[i-1]
         }
      
         if ((i%%100)==0)
          print(i/(burnin+thin*int))
     } 
    
     mucb=array(0,c(int));
     lamdacb=array(0,c(int));
     j=1
        
     for (i in (burnin+1):(burnin+thin*int)){
          muest=rbeta(1,(((mumc[i-1]^2)*(1-mumc[i-1]))/Vmu)+((mumc[i-1])/(1-mumc[i-1]))*(mumc[i-1]-1),
          ((((mumc[i-1])*((mumc[i-1])^2))/Vmu)-1+mumc[i-1]))   
          lambest=rgamma(1,lamdamc[i-1]^2/Vlamb,lamdamc[i-1]/Vlamb)
        
          alpha=exp(lpost(dados,muest,lambest)-lpost(dados,mumc[i-1],lamdamc[i-1]))
          alpha=alpha/(dbeta(muest,(((mumc[i-1]^2)*(1-mumc[i-1]))/Vmu)+((mumc[i-1])/(1-mumc[i-1]))*(mumc[i-1]-1),
          ((((mumc[i-1])*((mumc[i-1])^2))/Vmu)-1+mumc[i-1])))
          alpha=alpha/(dgamma(lambest,lamdamc[i-1]^2/Vlamb,lamdamc[i-1]/Vlamb))
          alpha=alpha*dbeta(mumc[i-1],((((muest^2)*(1-muest))/Vmu)+(muest/(1-muest))*(muest-1)),
          (((muest*(muest^2))/Vmu)-1+muest))
          alpha=alpha*dgamma(lamdamc[i-1],lambest^2/Vlamb,lambest/Vlamb)
        
          u=runif(1)
          if (u<alpha){
             mumc[i]=muest
             lamdamc[i]=lambest
          }
           else{
               mumc[i]=mumc[i-1]
               lamdamc[i]=lamdamc[i-1]
           }

          if ((i%%thin)==0){
              mucb[j]=mumc[i]
              lamdacb[j]=lamdamc[i]
              j=j+1
          }
          if ((i%%100)==0)
             print(i/(burnin+thin*int))
     }
     cbind(mucb,lamdacb)
}

x=rsimplex(1000,0.5,2)
result=posterior.simplex(x,1000)

