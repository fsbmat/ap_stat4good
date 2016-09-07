##############################################################################
###                                                                        ###
##           Geração de valores e estimação de parâmetros                   ##
###                            BS(alpha,beta)                              ###  
##                                                                          ##
##############################################################################

#geração de valores de uma distribuição Birnbaum-Saunders(alpha,beta)

#Observe que dentro da função log de verossimilhança foi usado a exp(lalpha) e exp(lbeta)
#onde lalpha=log(alpha) e lbeta=log(beta), isso foi feito para não existir problemas na otimização
#feita pela função optim que sem esses argumentos desta forma, utilizando alpha e beta, dá um erro
#pois as vezes considera log de números negativos. Ao fazer esta mudança não altera-se em nada a
#função pois a substituição usa exp(log(alpha))=alpha e exp(log(beta))=beta, isso é um macete
#para o R não acusar erro.
#Outra forma de resolver seria manter alpha e beta e acrescentar a função optim o comando 1, abaixo:
#Neste caso vou esconder os avisos. Eles continuarão lá, mas ocultos:
#1-suppressWarnings(optim(start,fn=Loglik,method="BFGS",hessian=T)$par)
#Outra forma é acrescentar anter do log de verossimilhança o código:
#if (any(c(t, beta, alpha, pi) < 0)) return(NA)
rm(list=ls())
cat("\014")

n<-100

#Gerando uma variável aleatória t com distribuição Birnbaum-Saunders(alpha,beta)
alpha=1
beta=2

#Note que teremos N amostras de tamanho 100 distintas, portanto t
#deve estar dentro do loop do Monte Carlo, alpha e beta deve estar fora do
#Loop pois os mesmos são fixos para todas as amostras.
N=1000
M=matrix(nrow=N,ncol=2)
for(i in 1:N){
  z<-rnorm(n,0,1)
  #Gerando uma variável aleatória t com distribuição Birnbaum-Saunders(alpha,beta)
  t<-cbind(beta*((alpha*z/2)+sqrt((alpha*z/2)^2+1))^2)
  #Verossimilhança
  Loglik<-function(par,dados){
    lalpha=par[1]
    lbeta=par[2]
    t<-dados
    ll<-sum(log(t+exp(lbeta)))-n*log(2*exp(lalpha))-(n/2)*(log(2*pi*exp(lbeta)))-(3/2)*sum(log(t))-((1/(2*exp(lalpha)^2))*sum((t/exp(lbeta))+(exp(lbeta)/t)-2))
    return(-ll)
  }
  
  #Utilizando a verossimilhança e a função optim para estimar os parâmetros alpha e beta que deram
  #origem as observações t observadas. Note que foi utilizado o log dos valores de chutes iniciais
  
  lalpha_0=log(2)
  lbeta_0=log(2)
  start=c(lalpha_0,lbeta_0)
  
  M[i,]=exp(optim(start,fn=Loglik,method="BFGS",dados=t,hessian=T)$par)
  #M[i,]=exp(nlminb(start,Loglik,dados=t)$par)
  
}
mean=colMeans(M)
mean

##############################################################################
###                                                                        ###
##                 Gráficos da função log de verossimilhança                ##
###                            BS(alpha,beta)                              ###  
##                                                                          ##
##############################################################################

#Usando o plot3d
require(rgl) # Para fazer os gráficos com rotação

alpha= seq(0.5, 3,l=n)
beta=seq(0.5,8,l=n)

# A função de verossimilhança para fazer o gráfico:
f<-function(alpha,beta){
  z1<-rnorm(n,0,1)
  t<-beta*((alpha*z1/2)+sqrt((alpha*z1/2)^2+1))^2
  (sum(log(t+beta))-n*log(2*alpha)-(n/2)*(log(2*pi*beta))-(3/2)*sum(log(t))-((1/(2*alpha^2))*sum((t/beta)+(beta/t)-2)))
}
z <- outer(alpha,beta, f)
#open3d()
#bg3d("white")
#material3d(col = "black")
persp3d(alpha, beta, z, col = "lightblue",
        xlab = "Shape", ylab = "Scale", zlab = "Loglik")

#agora vou plotar o gráfico da log de verossimilhança e visualizar o plot dos contornos/isovalores e
#o mapa de cores

#Neste caso posso considerar a log de verossimilhança sem o uso de artificíos para alpha e beta
#Log de Verossimilhança
LogLik<-function(par,dados){
  alpha=par[1]
  beta=par[2]
  ll<-sum(log(t+beta))-n*log(2*alpha)-(n/2)*(log(2*pi*beta))-(3/2)*sum(log(t))-((1/(2*alpha^2))*sum((t/beta)+(beta/t)-2))
  return(ll)
}

alpha=mean[1]
beta=mean[2]
pars.MV<-c(alpha,beta)

#criamos uma sequência adequada de pares de valores de (alpha,beta) e calculamos l(alpha,beta) 
#para cada um dos pares.

par.vals<-expand.grid(alpha= seq(0.5, 3,l=100),beta=seq(0.5,8,l=100)) 
dim(par.vals)
head(par.vals)

par.vals$logL<-apply(par.vals, 1, LogLik, dados = t) 
head(par.vals)
tail(par.vals)

#Note na sintaxe acima que a função apply aplica a função logveroN a cada par de valores 
#em cada linha de par.vals. Ao final o objeto |par.vals| contém na terceira coluna os valores 
#da log-verossimilhança correspondentes as valores dos parâmetros dados na primeira e segunda 
#colunas.

#superfície 3D gerada pela função persp()
with(par.vals, persp(unique(alpha), unique(beta), matrix(logL, 
                                                         ncol = length(unique(beta))),xlab = expression(alpha), 
                     ylab = expression(beta), zlab = expression(l(alpha,beta)), theta = 0, phi = 0,col = "lightblue")) 

# repeat{
#   for (i in 1:360) {
#     persp(x, y, z, theta = i, phi = 30, expand = 0.5, col = "lightblue")
#     }
# }

#help(persp)
#mapa de curvas de isovalores obtido com image()
with(par.vals,contour(unique(alpha),unique(beta), matrix(logL, ncol=length(unique(beta))),xlab = expression(alpha), ylab = expression(beta), nlev = 80)) 
points(pars.MV[1], pars.MV[2], pch = 4, cex = 1.5) 

#mapa de cores correspondentes aos valores gerado por image()
with(par.vals, image(unique(alpha), unique(beta), matrix(logL, ncol = length(unique(beta))),xlab = expression(alpha), ylab = expression(beta), col = gray(seq(0,1, length = 30)))) 
points(pars.MV[1], pars.MV[2], pch = 4, cex = 3)

#Para obtenção da função foi necessário especificar faixas de valores para alpha e beta. 
#A definição desta faixa foi feita após várias tentativas pois depende do problema, em especial
#do número e variabilidade dos dados.
#as funções gráficas utilizadas requerem: dois vetores de tamanhos n1 e n2 com os valores dos 
#argumentos da função e os valores da função em uma matrix de dimensão n1 × n2. Por isto usamos 
#unique() para extrair os valores dos argumentos, sem repeti-los e matrix() para os valores da 
#função.
#na função perp() as argumentos theta e phi são utilizados para rotacionar o gráfico a fim de 
#se obter uma melhor visualização.
#o valor das estimativas de máxima verossimilhança são indicados por x nos dois últimos gráficos. 
#Neste caso eles foram fixados no objeto pars.MV, mas eles podem ser obtidos analiticamente. 
