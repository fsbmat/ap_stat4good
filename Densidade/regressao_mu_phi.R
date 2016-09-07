##############################################################################
###                                                                        ###
##    Criando uma estrutura de regressão do tipo beta0+beta1*X para mu      ##
###                     e alpha0+alpha1*Z para phi                          ##
##                             BS(mu,phi)                                  ###  
##                                                                          ##
##############################################################################
rm(list=ls())
cat("\014")

N=100
#m e m1 são as matrizes que receberão as estimativas no final do processo de 
#estimação
m=matrix(ncol=4,nrow=N)
m1=matrix(ncol=4,nrow=N)
#Valores iniciais dos parâmetros usados para gerar t
#Ou seja, Valor verdadeiro dos parâmetros
beta0=2
beta1=-1
alpha0=3
alpha1=1
truevalue=c(beta0,beta1,alpha0,alpha1)
#Tamanho das amostras
n=100
#Vetor de parâmetros
beta=matrix(c(beta0,beta1),nrow=2,ncol=1)
alpha=matrix(c(alpha0,alpha1),nrow=2,ncol=1)
#Vetor de 1's
const1 <- rep(1,n);
const <- cbind(const1);
#Vetor de cováriavel com distribuição unif(0,1)
X1=matrix(runif(n, 0, 1),nrow=n,ncol=1)
Z1=matrix(runif(n, 0, 1),nrow=n,ncol=1)
#Matrizes de covariáveis
X <-matrix(c(const,X1),nrow=n,ncol=ncol(X1)+1)
Z <-matrix(c(const,Z1),nrow=n,ncol=ncol(Z1)+1)
#Número de colunas de X e Z
p=ncol(X)
q=ncol(Z)
#Vetor de médias
mu=exp(X%*%beta)
#Vetor de dispersão
phi=exp(Z%*%alpha)

#Gerando uma variável aleatória t com distribuição Birnbaum-Saunders(mu,phi)
remove(beta,beta0,beta1,alpha,alpha0,alpha1)
#Monte Carlo para estimação dos parâmetros beta0, beta1 e beta2

for (i in 1:N){
  #set.seed(123)  
  z<-cbind(rnorm(n,0,1))
  t<-((phi*mu)/(phi+1))*((z/sqrt(2*phi))+(sqrt((z/sqrt(2*phi))^2+1)))^2
  
  
  #Um motivo de erro comum na função abaixo é esquecer o parentêses dentro do colchete
  #de theta[(p+1):(p+q)]  
  Loglik<-function(theta,dados){
    p=ncol(X)
    q=ncol(Z)
    beta=theta[1:p]
    alpha=theta[(p+1):(p+q)]
    mu=exp(X%*%beta)
    phi=exp(Z%*%alpha)
    lv=sum((phi/2)-(3/2)*log(t)+(1/2)*log(phi+1)-log(4*sqrt(pi*mu))+log(t+(phi*mu/(phi+1)))-(phi/4)*((t*(phi+1)/(phi*mu))+(phi*mu/(t*(phi+1)))))
    #lv=sum(log(((exp(phi/2)*sqrt(phi+1))/(4*sqrt(pi*mu)*t^(3/2)))*(t+((phi*mu)/(phi+1)))*exp((-phi/4)*((t*(phi+1)/(phi*mu))+(phi*mu/(t*(phi+1)))))))
    return(-lv)
  }
  #Chute inicial para as funções de estimação
  start=cbind(5,3,1,-1)
  
  #alpha=c(1,2)
  #beta=c(2,1)
  
  #Estimation with function optim
  bs_op=optim(start,Loglik,method="BFGS",dados=t,hessian = T) 
  m[i,]=bs_op$par
  bs_nl=nlminb(start, Loglik,dados=t)
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
         eqm(x=m[,2],bs_op=truevalue[2]),
         eqm(x=m[,3],bs_op=truevalue[3]),
         eqm(x=m[,4],bs_op=truevalue[4]))

#Estimated mean squared error of each parameter
eqmest1=c(eqm1(x=m1[,1],bs_nl=truevalue[1]),
          eqm1(x=m1[,2],bs_nl=truevalue[2]),
          eqm1(x=m1[,3],bs_nl=truevalue[3]),
          eqm1(x=m1[,4],bs_nl=truevalue[4]))


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
hist(m[,3],prob=T);
rug(m[,3])
curve(expr = dnorm(x,mean=mean(m[,3]),sd=sd(m[,3])),add=T, col="red")
hist(m[,4],prob=T);
rug(m[,4])
curve(expr = dnorm(x,mean=mean(m[,4]),sd=sd(m[,4])),add=T, col="red")

hist(m1[,1],prob=T);
rug(m1[,1])
curve(expr = dnorm(x,mean=mean(m1[,1]),sd=sd(m1[,1])),add=T, col="red")
hist(m1[,2],prob=T);
rug(m1[,2])
curve(expr = dnorm(x,mean=mean(m1[,2]),sd=sd(m1[,2])),add=T, col="red")
hist(m1[,3],prob=T);
rug(m1[,3])
curve(expr = dnorm(x,mean=mean(m1[,3]),sd=sd(m1[,3])),add=T, col="red")
hist(m1[,4],prob=T);
rug(m1[,4])
curve(expr = dnorm(x,mean=mean(m1[,4]),sd=sd(m1[,4])),add=T, col="red")
