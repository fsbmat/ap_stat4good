###########################################################################
###                                                                     ###
##                 Gráficos da Densidade BS(mu,phi)                      ##
###                      nova parametrização                            ###
###########################################################################

rm(list=ls())
cat("\014")

#Para criar uma estrutura de regressão é interessante que se faça uma reparametrização
#de forma que E(T)=mu, para que façamos a modelagem da média da distribuição. Assim,
#considere:

#alpha=sqrt(2/phi) e beta=(phi*mu)/(phi+1) \rightarrow phi=(2/alpha^2) e mu=beta*(1+((alpha^2)/2))
#phi>0 e mu>0

#A nova densidade é:

bs<-function(t,mu,phi){
  fdp=((exp(phi/2)*sqrt(phi+1))/(4*sqrt(pi*mu)*t^(3/2)))*(t+((phi*mu)/(phi+1)))*exp((-phi/4)*((t*(phi+1)/(phi*mu))+(phi*mu/(t*(phi+1)))))
  return(fdp)
}

#Teste para ver se a integral da densidade é igual a 1.
mu=1
phi=2
integrand <- function(t){((exp(phi/2)*sqrt(phi+1))/(4*sqrt(pi*mu)*t^(3/2)))*(t+((phi*mu)/(phi+1)))*exp((-phi/4)*((t*(phi+1)/(phi*mu))+(phi*mu/(t*(phi+1)))))}
integrate(integrand, lower = 0, upper = Inf)

#Note que com a mudança de mu, mantendo phi fixo, há uma alteração na média e na variança da variável aleatória
#Com o aumento de beta há um aumento da média e da variança
t=seq(0.5,6,by=0.01)
plot(t,bs(t,1,100),main='T ~ Birnbaum-Saunders(mu, phi=1)',ylab='Densidade',type='l')
lines(t,bs(t,1.5,100),col=2,lty=2, lwd=1)
lines(t,bs(t,2,100),col=3,lty=3, lwd=2)
lines(t,bs(t,2.5,100),col=4,lty=4, lwd=3)
lines(t,bs(t,3,100),col=5,lty=5, lwd=3)
lines(t,bs(t,3.5,100),col=6,lty=6, lwd=2)
legend(4.3, 2.7, c("BS(t,1,100)","BS(t,1.5,100)","BS(t,2,100)","BS(t,2.5,100)","BS(t,3,100)","BS(t,3.5,1)"), fill=1:6)

#Enquanto que ao manter mu fixo e mudar o phi há uma mudança na assimetria do gráfico
t=seq(0,4,by=0.01)
plot(t,bs(t,1,2),ylim = c(0,2.8),main='T ~ Birnbaum-Saunders(mu=1, phi)',ylab='Densidade',type='l')
lines(t,bs(t,1,5),col=2,lty=2, lwd=1)
lines(t,bs(t,1,10),col=3,lty=3, lwd=2)
lines(t,bs(t,1,25),col=4,lty=4, lwd=3)
lines(t,bs(t,1,50),col=5,lty=5, lwd=3)
lines(t,bs(t,1,100),col=6,lty=6, lwd=2)
legend(2.5, 2.5, c("BS(t,1,2)","BS(t,1,5)","BS(t,1,10)","BS(t,1,25)","BS(t,1,50)","BS(t,1,100)"), fill=1:6)

V=function(mu,phi){
  vt=2*(mu^2)*(2*phi+5)/((phi+1)^2)
  return(vt)
}

#Ao manter mu fixo e aumentar phi a variância de T (Var(T)) tende a zero 
phi=seq(0.001,100,by=0.01)
plot(phi,V(0.6,phi),ylim = c(0,2.8),main='T ~ Birnbaum-Saunders(mu=1, phi)',ylab='Var(T)',type='l')
lines(phi,V(1,phi),col=2,lty=2, lwd=1)
lines(phi,V(1.5,phi),col=3,lty=3, lwd=2)
lines(phi,V(2,phi),col=4,lty=4, lwd=3)
lines(phi,V(2.5,phi),col=5,lty=5, lwd=3)
lines(phi,V(3,phi),col=6,lty=6, lwd=2)
legend(60, 2.5, c("Var(T;0.6,phi)","Var(T;1,phi)","Var(T;1.5,phi)","Var(T;2,phi)","Var(T;2.5,phi)","Var(T;3,phi)"), fill=1:6)

