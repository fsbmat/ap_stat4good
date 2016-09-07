rm(list=ls())
cat("\014")
par(mfrow=c(1,1))

###########################################################################
###                                                                     ###
##                 Gráficos da Densidade BS(alpha,beta)                  ##
###                     parametrização original                         ###
###########################################################################

#Função distribuição acumulada da Birnbaum-Saunders com parâmetros de forma, alpha, e de escala beta.
#Shape=alpha, scala=beta

acbs<-function(t,alpha,beta){
  ff<-pnorm((1/alpha)*(sqrt(t/beta)-sqrt(beta/t)))
  return(ff)
}

#Note que com a mudança de alpha, mantendo beta fixo, há uma alteração na assimetria do gráfico
t=seq(0,3,by=0.01)
plot(t,acbs(t,0.1,1),main='T ~ Birnbaum-Saunders(shape, scale=1)',ylab='F(T<=t)',type='l')
lines(t,acbs(t,0.3,1),col=2,lty=2, lwd=1)
lines(t,acbs(t,0.5,1),col=3,lty=3, lwd=2)
lines(t,acbs(t,0.75,1),col=4,lty=4, lwd=3)
lines(t,acbs(t,1,1),col=5,lty=5, lwd=3)
lines(t,acbs(t,1.5,1),col=6,lty=6, lwd=2)
legend(2.3, 0.55, c("BS(t,0.1,1)","BS(t,0.3,1)","BS(t,0.5,1)","BS(t,0.75,1)","BS(t,1,1)","BS(t,1.5,1)"), fill=1:6)

#Função densidade da Birnbaum-Saunders com parâmetros de forma, alpha, e de escala beta.
#Shape=alpha, scala=beta

bs<-function(t,alpha,beta){
  f<-((t+beta)/(2*alpha*sqrt(2*pi*beta)))*(t^(-3/2))*exp((-1/(2*alpha^2))*((t/beta)+(beta/t)-2))
  return(f)
}

alpha=10
beta=2
integrand <- function(t){((t+beta)/(2*alpha*sqrt(2*pi*beta)))*(t^(-3/2))*exp((-1/(2*alpha^2))*((t/beta)+(beta/t)-2))}
integrate(integrand, lower = 0, upper = Inf)

#Note que com a mudança de alpha, mantendo beta fixo, há uma alteração na assimetria do gráfico
t=seq(0,3,by=0.01)
plot(t,bs(t,0.1,1),main='T ~ Birnbaum-Saunders(shape, scale=1)',ylab='Densidade',type='l')
lines(t,bs(t,0.3,1),col=2,lty=2, lwd=1)
lines(t,bs(t,0.5,1),col=3,lty=3, lwd=2)
lines(t,bs(t,0.75,1),col=4,lty=4, lwd=3)
lines(t,bs(t,1,1),col=5,lty=5, lwd=3)
lines(t,bs(t,1.5,1),col=6,lty=6, lwd=2)
legend(2.3, 4, c("BS(t,0.1,1)","BS(t,0.3,1)","BS(t,0.5,1)","BS(t,0.75,1)","BS(t,1,1)","BS(t,1.5,1)"), fill=1:6)

#Enquanto que ao manter alpha fixo e mudar o beta há uma mudança na média e na variança da variável aleatória
#Com o aumento de beta há um aumento da média e da variança
t=seq(0.2,3.5,by=0.01)
plot(t,acbs(t,0.1,0.5),main='T ~ Birnbaum-Saunders(shape=0.1, scale)',ylab='F(T<=t)',type='l')
lines(t,acbs(t,0.1,0.75),col=2,lty=2, lwd=1)
lines(t,acbs(t,0.1,1),col=3,lty=3, lwd=2)
lines(t,acbs(t,0.1,1.5),col=4,lty=4, lwd=3)
lines(t,acbs(t,0.1,2),col=5,lty=5, lwd=3)
lines(t,acbs(t,0.1,2.5),col=6,lty=6, lwd=2)
legend(2.55, 0.55, c("BS(t,0.1,0.5)","BS(t,0.1,0.75)","BS(t,0.1,1)","BS(t,0.1,1.5)","BS(t,0.1,2)","BS(t,0.1,2.5)"), fill=1:6)

#Enquanto que ao manter alpha fixo e mudar o beta há uma mudança na média e na variança da variável aleatória
#Com o aumento de beta há um aumento da média e da variança
t=seq(0,3,by=0.01)
plot(t,bs(t,0.1,0.5),main='T ~ Birnbaum-Saunders(shape=0.1, scale)',ylab='Densidade',type='l')
lines(t,bs(t,0.1,0.75),col=2,lty=2, lwd=1)
lines(t,bs(t,0.1,1),col=3,lty=3, lwd=2)
lines(t,bs(t,0.1,1.5),col=4,lty=4, lwd=3)
lines(t,bs(t,0.1,2),col=5,lty=5, lwd=3)
lines(t,bs(t,0.1,2.5),col=6,lty=6, lwd=2)
legend(2.2, 8, c("BS(t,0.1,0.5)","BS(t,0.1,0.75)","BS(t,0.1,1)","BS(t,0.1,1.5)","BS(t,0.1,2)","BS(t,0.1,2.5)"), fill=1:6)
