
rm(list=ls(all=TRUE))
set.seed(7)

# Função que gera valores da Burr XII com dois parâmetros
rburrxii=function(n,p,b){
  u=runif(n)
  p=5
  b=5
  x=(((1-u)^(-1/p))-1)^(1/b)
  return(x)
}

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

p=4;b=5

x=rburrxii(25,p,b)
result=posterior.burr(x,1000,theta=c(p,b))
