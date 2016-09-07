##############################################################################
###                                                                        ###
##    Criando uma estrutura de regressão do tipo beta0+beta1*X para mu      ##
###                            BS(mu,phi)                                  ###  
##                                                                          ##
##############################################################################
rm(list=ls())
cat("\014")

N=1000
#m,m1,m2 e m3 são as matrizes que receberão as estimativas no final do processo de 
#estimação
m=matrix(ncol=2,nrow=N)
m1=matrix(ncol=2,nrow=N)
#Valores iniciais dos parâmetros usados para gerar t
#Ou seja, Valor verdadeiro dos parâmetros
beta0=2
beta1=1
truevalue=c(beta0,beta1)
#Tamanho das amostras
n=100
#Vetor de parâmetros
beta=matrix(c(beta0,beta1),nrow=2,ncol=1)
#Vetor de 1's
const1 <- rep(1,n);
const <- cbind(const1);
#Vetor de cováriavel com distribuição unif(0,1)
X1=matrix(runif(n, 0, 1),nrow=n,ncol=1)
#Matriz de covariáveis
X <-matrix(c(const,X1),nrow=n,ncol=ncol(X1)+1)
#Número de colunas de X
p=ncol(X)
#Vetor de médias
mu=exp(X%*%beta)

#Gerando uma variável aleatória t com distribuição Birnbaum-Saunders(mu,phi)
phi=2
remove(beta,beta0,beta1)
#Monte Carlo para estimação dos parâmetros beta0, beta1 e beta2

for (i in 1:N){
  #set.seed(123)  
  z<-cbind(rnorm(n,0,1))
  t<-((phi*mu)/(phi+1))*((z/sqrt(2*phi))+(sqrt((z/sqrt(2*phi))^2+1)))^2
  
  Loglik<-function(beta,dados){
    p=ncol(X)
    mu=exp(X%*%beta[1:p])
    phi=phi
    lv=sum((phi/2)-(3/2)*log(t)+(1/2)*log(phi+1)-log(4*sqrt(pi*mu))+log(t+(phi*mu/(phi+1)))-(phi/4)*((t*(phi+1)/(phi*mu))+(phi*mu/(t*(phi+1)))))
    #lv=sum(log(((exp(phi/2)*sqrt(phi+1))/(4*sqrt(pi*mu)*t^(3/2)))*(t+((phi*mu)/(phi+1)))*exp((-phi/4)*((t*(phi+1)/(phi*mu))+(phi*mu/(t*(phi+1)))))))
    return(-lv)
  }
  #Chute inicial para as funções de estimação
  start=c(10,20)
  
  #Estimation with function optim
  bs_op=optim(start,Loglik,method="BFGS",dados=t,hessian = T) 
  m[i,]=bs_op$par
  bs_nl=nlminb(start, Loglik)
  m1[i,]=bs_nl$par
}
mest=colMeans(m)
mest1=colMeans(m1)

#calculating the standard deviation of each column of the array of parameters m
dest=apply(m,2,sd)
dest1=apply(m1,2,sd)

#root mean square error in the calculation of each column of the array of parameters m in relation to the true value of the parameter
eqm=function(x,bs_op){ 
  k=length(x)
  sqrt(sum(((x-bs_op)^2))/k)}

eqm1=function(x,bs_nl){ 
  k=length(x)
  sqrt(sum(((x-bs_nl)^2))/k)}

#Estimated mean squared error of each parameter 

eqmest=c(eqm(x=m[,1],bs_op=truevalue[1]),
         eqm(x=m[,2],bs_op=truevalue[2]))

#Estimated mean squared error of each parameter
eqmest1=c(eqm1(x=m1[,1],bs_nl=truevalue[1]),
          eqm1(x=m1[,2],bs_nl=truevalue[2]))


# Table with the true values of the parameters and the average
# Standard deviation and mean square error of the estimated parameters
tab=data.frame(truevalue,mean=mest,sd=dest,eqm=eqmest)
tab1=data.frame(truevalue,mean=mest1,sd=dest1,eqm=eqmest1)

tab
tab1

par(mfrow=c(2,2))
hist(m[,1],prob=T);
rug(m[,1])
curve(expr = dnorm(x,mean=mean(m[,1]),sd=sd(m[,1])),add=T, col="red")
hist(m[,2],prob=T);
rug(m[,2])
curve(expr = dnorm(x,mean=mean(m[,2]),sd=sd(m[,2])),add=T, col="red")
hist(m1[,1],prob=T);
rug(m1[,1])
curve(expr = dnorm(x,mean=mean(m1[,1]),sd=sd(m1[,1])),add=T, col="red")
hist(m1[,2],prob=T);
rug(m1[,2])
curve(expr = dnorm(x,mean=mean(m1[,2]),sd=sd(m1[,2])),add=T, col="red")
