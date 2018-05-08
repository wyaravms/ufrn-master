
# exercicio inferencia

# 1) Fixe  \theta_0 para algum f(x|\theta)
theta0 = 2

# 2) criar sequencia de valores para \theta
theta = seq(1,5, 0.01)

# 3) fixar n e gerar X de X~f(x|\theta)
n = 10
x = rexp(n, theta0)

# obter L(\theta, x) para da \theta seq
lexp - function(x, n, theta){(theta^n)*exp(-theta*sum(x))}
loglexp = function(x, n, theta){(n*log(theta))-(theta*sum(x))}
l=loglexp(x, n, theta)

# representar graficamente \theta x L(\theta, x) para \theta seq e comentar
#par(mfrow=c(1,2))
plot(theta, l, type="l", col=2)

# function geral para gráfico da log verossimilhança
par(mfrow=c(2,3))
for (i in 1:6){
  theta0=2
  theta=seq(1,5, 0.01)
  n=10^i
  x=rexp(n, theta0)
  loglexp<-function(x, n, theta){(n*log(theta))-(theta*sum(x))}
  l=loglexp(x, n, theta)
  plot(theta, l, type="l", col=2, main=paste("n= ", n))
}

# replicas - primeira tentantiva 
replicas = function(n){
nrep = 1000
theta0 = 2
x = NULL
a = 0
b = 0
for (i in 1:nrep){
 theta = seq(1,10, 0.1)
 x = rexp(n,theta0)
 #lexp = function(x, n, theta){(n*log(theta))-(theta*sum(x))}
 loglexp = function(x, n, theta){(n*log(theta))-(theta*sum(x))}
 l = loglexp(x, n, theta0)
 l1 = loglexp(x, n, sample(theta,1, replace = FALSE))
 if (l>l1)
  {a  = a+1} 
    else {b = b+1}
}
prob = (1/nrep)*sum(a)
prob
}

replicas(10)
replicas(100)
replicas(1000)
replicas(10000)
