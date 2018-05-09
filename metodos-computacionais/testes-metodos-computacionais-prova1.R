
# ---------------------------------------------------------------#
# Geradores Congruenciais
# ---------------------------------------------------------------#

# encontrar a que forneca periodo completo
b=NULL
for(a in 2:10){
      b[a]=(a^10-1)%%11
      a
}

# ---------------------------------------------------------------#
# Método da tranformada inversa para "Distribuicoes continuas"
# ---------------------------------------------------------------#

# generate random values
# exponential distribution
n=1000
lambda=2
     
U<-runif(n)
x<-(-log(1-U)/lambda)
x

#outra parametrizacao
x2<-(-lambda*log(1-U))
x2      
  
y<-rexp(n,lambda)

par(mfrow=c(1,3))
hist(x)
hist(x2)
hist(y)

#curva da densidade
d=density(x)
d2=density(x2)
dy=density(y)
plot(d, col="blue")
lines(d2, col="red")
lines(dy, col="yellow") 

# pareto distribution
n=1000
beta=2
k=2    
     
U<-runif(n)
x<-(beta/((1-U)^(1/k)))
x     

# package https://cran.r-project.org/web/packages/VGAM/ 
library(VGAM)
y<-rpareto(n,beta,k)

par(mfrow=c(1,2))
hist(x)
hist(y)

# burr XII distribution com parametro de escala igual a 0.
n=100
c=1
k=1
alpha=1

U<-runif(n)
x<-((alpha)*(((1/(1-U)^(1/k))-1)^(1/c)))
x 

fburrxii<-function(x, c, k, alpha)
{
d<-(((k*c)/alpha)*((x/alpha)^(c-1)))/((1+(x/alpha)^c)^(k+1))
#d<-k*c*x^(c-1)*(1+x^c)^(-k-1)
d
}

hist(x, prob=TRUE)
curve(fburrxii(x, c, k, alpha), add=TRUE, col=2)
library(VaRES)
curve(dburr7(x, k, c, log=FALSE), add=TRUE, col=4)

# ---------------------------------------------------------------#
# metodo da transformada inversa para "variavel discreta"
# ---------------------------------------------------------------#

amostra<-NULL
for(k in 1:100){
u<-runif(1,0,1)
if(u<=0.2) amostra[k]<-1
if((0.2<u)&(u<=0.35))amostra[k]<-2
if((0.35<u)&(u<=0.6))amostra[k]<-3
if(u>0.6)amostra[k]<-4
}

# using acumulate
amostra<-NULL
for(k in 1:100){
u<-runif(1,0,1)
if(u<0.2) amostra[k]<-1
if(u<0.35)amostra[k]<-2
if(u<0.6)amostra[k]<-3
if(u<1)amostra[k]<-4
}

# using floor to generate uniform distribution discrete
n<-10
u<-runif(100,0,1)
amostra<-floor(n*u)+1


# ---------------------------------------------------------------#
# metodo da transformada inversa para "variavel discreta"
# ---------------------------------------------------------------#

# generate randon permutation
n=10
p=seq(1,n,1)
k=n
while(k>1){
u<-runif(1,0,1)
I<-floor(n*u)+1
x<-p[I]
y<-p[k]
p[I]<-y
p[k]<-x
k=k-1
}

#gerar normal padrao
x<-NULL
z<-NULL
for(i in 1:100){
u1<-runif(1,0,1)
y<-(-log(u1))
u2<-runif(1,0,1)
if (u2<=(exp(-((y-1)^2/2)))) {x[i]=y
u3<-runif(1,0,1) 
if (u3<=0.5) z[i]=x[i] else z[i]=-x[i]
}
}

# outro maneira 
x<-NULL
z<-NULL
for(i in 1:100){
u1<-runif(1,0,1)
y1<-(-log(u1))
u2<-runif(1,0,1)
y2<-(-log(u2))
if (y2>=(((y1-1)^2/2))) {x[i]=y1
u3<-runif(1,0,1) 
if (u3<=0.5) z[i]=x[i] else z[i]=-x[i]
}
}

# metodo inversa discreta
amostra<-NULL
for(k in 1:100){
u<-runif(1,0,1)
if(u<=0.2)amostra[k]<-1
if((0.2<u)&(u<=0.35)) amostra[k]<-2
if((0.35<u)&(u<=0.6)) amostra[k]<-3
if(0.6<u) amostra[k]<-4
}
amostra

# outra maneira inversa discreta FAIL
amostra<-NULL
for(k in 1:100){
u<-runif(1,0,1)
if(u<0.2)amostra[k]<-1  
break
if(u<0.35) amostra[k]<-2    
break
if(u<0.6) amostra[k]<-3
break
if(0.6<u) amostra[k]<-4
}
amostra


# gerar ditribuicao discreta Y pelo metodo da rejeicao
n<-100
p<-c(0.11,0.12,0.09,0.08,0.12,0.10,0.09,0.09,0.10,0.10)
amostra<-NULL
k<-0
while(k<=n-1){
u1<-runif(1,0,1)
x<-floor(10*u1)+1
pass<-10*p[x]/1.2
if(u1<=pass)((k<-k+1)&(amostra[k]<-x))
}

# geracao de ocorrencia da Distribuicao Poisson
Rpois<-function(n,lambda){
x=NULL    

for(j in 1:n){ 
p=exp(-lambda)
i=0
F=p
u<-runif(1,0,1)                                              
      while (u>F) { 
            p=(lambda*p)/(i+1)
            F=F+p
            i=i+1
      }
      x[j]=i         
}
return(x)
}

# geracao de ocorrencia da Distribuicao Binomial
Rbinom<-function(n,size,p){
x=NULL

for (j in 1:n){
c=p/(1-p)
i=0
pr=(1-p)^size
F=pr
u<-runif(1,0,1)
    while(u>F){
           pr=((c*(size-i)/(i+1)))*pr
           F=F+pr
           i=i+1
    }
    x[j]=i         
}
return(x)
}

# ---------------------------------------------------------------#
# metodo da composicao
# ---------------------------------------------------------------#

for (i in 1:1000){
u1<-runif(1,0,1)
u2<-runif(1,0,1)
if (u1<0.5) {x[i]=floor(10*u1)+1} else {x[i]=floor(5*u2)+6}
}

# ---------------------------------------------------------------#
# gerador variaveis da gamma distribution
# ---------------------------------------------------------------#

Rgamma<-function(m, n, lambda){
u=NULL
y=NULL

for(j in 1:m){
for(i in 1:n){
u[i]<-runif(1,0,1)
}
y[j]<-(-lambda*log(prod(u)))
}
return(y)
}

n=10
m=100
lambda=2
y=Rgamma(m, n, lambda)

d=density(y)
x=rgamma(100,n,1/lambda)
d1=density(x)
plot(d, col="red")
lines(d1, col="blue")

# ---------------------------------------------------------------#
# gerador variaveis da standard normal distribution  Box-Miller
# ---------------------------------------------------------------#
 
RnormBM<-function(n){
x=NULL
y=NULL

for (i in 1:n){
u1<-runif(1,0,1)
u2<-runif(1,0,1)
R2<-(-2*log(u1))
Theta<-2*pi*u2
x[i]<-(sqrt(R2)*cos(Theta))
y[i]<-(sqrt(R2)*sin(Theta))
}
cbind (x,y)
}

n=100
x=RnormBM(n)

d=density(x[,1])
d1=density(x[,2])

plot(d, col="red")
lines(d1, col="blue")

# ---------------------------------------------------------------#
# gerador variaveis da standard normal distribution  Polar-Method
# ---------------------------------------------------------------#

RnormPolar<-function(n){
x=NULL
y=NULL
i=0
while (i<n){
u1<-runif(1,0,1)
u2<-runif(1,0,1)
v1<-2*u1-1
v2<-2*u2-1
s=v1^2+v2^2
if(s<1) {          
   x[i]=sqrt((-2*log(s))/s)*v1
   y[i]=sqrt((-2*log(s))/s)*v2
   i=i+1
   }
}
cbind (x,y)
}

n=10000
x=RnormPolar(n)

d=density(x[,1])
d1=density(x[,2])

plot(d, col="red")
lines(d1, col="blue")



# ---------------------------------------------------------------#
# Metodos de otimizaçao
# ---------------------------------------------------------------#

#Simulate Annealing
x<-NULL
y<-NULL
max1<-7000 #Lk
constante<-0.4 #a
par.atual<-c(2,3) #theta_0

x[1]<-par.atual[1]
y[1]<-par.atual[2]
 
for (n in 1:max1){
    temp<-(constante)/(1000*log(n+1)) #Ck
    a<-par.atual[1]+runif(1,-0.05,0.05)
    b<-par.atual[2]+runif(1,-0.05,0.05)
    par.cand<-c(a,b)
    f.atual<- -exp(-(par.atual[1]^2 + par.atual[2]^2))
    f.cand <- -exp(-(par.atual[1]^2 + par.atual[2]^2))
    diferenca<- -(f.cand-f.atual)
    razao<- exp(diferenca/temp)
    fanta<- min(razao, 1)
    u<- runif(1,0,1)
    if(u<fanta) par.atual<-par.cand
    x[n]<-par.atual[1]
    y[n]<-par.atual[2]
    }

minimo<-c(x[7000], y[7000])
minimo   

# Exemplo para Inferencia 
#gerar os dados
n<-20
theta<-0.4
u<-runif(n,0,1)
dados<-(-1+2*(1/4-theta*(1/2-theta/4-u))^(0.5))/theta
#funcao escore e log verossimilhanca
e<-function(x){-n*log(2)+sum(log(1+x*dados))}
s<-function(y){sum((dados)/(1+y*dados))}
ds<-function(z){-sum(((dados)/(1+z*dados))^2)}

# Metodo Newton-Raphson
theta.0<-0.15
precisao<-0.000001
dif<-1
while(dif>precisao){
num<-s(theta.0)
den<-ds(theta.0)
theta.1<-theta.0-(num/den)
dif<-abs(theta.1-theta.0)
theta.0<-theta.1
}

raiz<-theta.0
raiz

# Metodo Escore Fisher
theta.0<-0.15
precisao<-0.00001
dif<-1
while(dif>precisao){
num<-s(theta.0)
a<-2*theta.0
b<-log((1+theta.0)/(1-theta.0))-a
den<-n*(1/(2*theta.0^3))*b
theta.1<-theta.0+(num/den)
dif<-abs(theta.1-theta.0)
theta.0<-theta.1
}

raiz<-theta.0
raiz

# optimize beta usando uniformes

alpha=2
beta=2
x=seq(0,1,0.01)
f<-function(x, alpha, beta){-((gamma(alpha+beta))/(gamma(alpha)*gamma(beta))*x^(alpha-1)*(1-x)^(beta-1))}
optimize(f,c(0,1), tol=0.0001, alpha=2, beta=2)
plot(f)





