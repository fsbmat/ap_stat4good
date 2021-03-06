Introdução
==========

<p style="text-align: justify;">
Neste post é utilizada a distribuição Birnbaum-Saunders para ilustrar
alguns conceitos de métodos computacionais em Inferência Estatística.
Essa distribuição foi proposta por Birnbaum e Saunders (1969), no artigo
intitulado **A new family of life distributions** e é, comumente,
considerado na modelagem de tempo de vida de materiais e equipamentos
sujeitos a cargas dinâmicas através de modelos de dano acumulado. Sendo,
amplamente, utilizada na área de engenharia, na indústria, em negócios,
na análise de confiabilidade, na análise de sobrevivência, em ciências
ambientais e ciências médicas e em diversas outras áreas, pois possui
propriedades interessantes e uma relação próxima com a distribuição
normal, o que a torna, do ponto de vista de aplicação, uma alternativa
mais atraente para as bem conhecidas distribuições Weibull,
log-logística, log-normal, gama e modelos inversos Gaussianos.
</p>
<p style="text-align: justify;">
Além do contexto inferencial, é discutido algumas observações sobre
questões numéricas quando trabalha-se com a função de verossimilhança na
perspectiva computacional fazendo uso do software `R`.
</p>
Função de Distribuição Acumulada Birnbaum-Saunders
--------------------------------------------------

<p style="text-align: justify;">
Seja *T* uma variável aleatória representando o tempo até a ocorrência
do evento de interesse, então assumindo que essa variável segue a
distribuição Birnbaum-Saunders, tem-se que a sua função de distribuição
acumulada (f.d.a.) é dada por:
</p>
\begin{equation}\label{fun1}
F_{T}(t;\alpha,\beta)=P(T\leq t)=\Phi\left[\dfrac{1}{\alpha}\left(\sqrt{\frac{t}{\beta}}-\sqrt{\frac{\beta}{t}}\right)\right],\ t>0
\end{equation}
<p style="text-align: justify;">
em que *Φ*(.) é a fda de uma distribuição normal padrão. Dizemos que *T*
segue uma distribuição BS, com parâmetros de forma *α* &gt; 0 e de
escala *β* &gt; 0, que é usualmente denotada por *T* ∼ *B**S*(*α*, *β*).
</p>
Função Densidade
----------------

<p style="text-align: justify;">
Considerando a distribuição acumulada da variável aleatória *T* dada em
$(\\ref{fun1})$ a sua correspondente função densidade de probabilidade
(fdp) é dada por
\begin{equation}\label{fun2}
f_{T}(t)=\dfrac{t^{-\frac{3}{2}}(t+\beta)}{2\sqrt{2\pi}\alpha\sqrt{\beta}}\exp\left[-\frac{1}{2\alpha^{2}}\left(\frac{t}{\beta}+\frac{\beta}{t}-2\right)\right]
\end{equation}
em que *α* &gt; 0 e *β* &gt; 0.
</p>
Função de distribuição acumulada, função densidade e seus respectivos gráficos no R.
------------------------------------------------------------------------------------

    #Shape=alpha, scala=beta
    acbs<-function(t,alpha,beta){
      ff<-pnorm((1/alpha)*(sqrt(t/beta)-sqrt(beta/t)))
      return(ff)
    }

    #Função densidade 
    bs<-function(t,alpha,beta){
      f<-((t+beta)/(2*alpha*sqrt(2*pi*beta)))*(t^(-3/2))*exp((-1/(
        2*alpha^2))*((t/beta)+
                       (beta/t)-2))
      return(f)
    }

<p style="text-align: justify;">
Conta simples para verificar se a integral da função densidade é igua a
1. Isso é somente uma forma de analisar se a digitação da densidade esta
correta.
</p>
    alpha=10
    beta=2
    integrand <- function(t){((t+beta)/(2*alpha*sqrt(2*pi*beta)))*(t^(-3/2))*
        exp((-1/(2*alpha^2))*((t/beta)+(beta/t)-2))}
    integrate(integrand, lower = 0, upper = Inf)

    ## 1 with absolute error < 4,1e-05

<p style="text-align: justify;">
Note que com a mudança de alpha, mantendo beta fixo, há uma alteração na
assimetria do gráfico, veja os gráficos da distribuição acumulada e da
densidade:
</p>
    t=seq(0,3,by=0.01)
    plot(t,acbs(t,0.1,1),main=expression(T%~%BS(alpha,beta==1)),ylab=expression(P(T<=t)==F(t)),type='l')
    lines(t,acbs(t,0.3,1),col=2,lty=2, lwd=1)
    lines(t,acbs(t,0.5,1),col=3,lty=3, lwd=2)
    lines(t,acbs(t,0.75,1),col=4,lty=4, lwd=3)
    lines(t,acbs(t,1,1),col=5,lty=5, lwd=3)
    lines(t,acbs(t,1.5,1),col=6,lty=6, lwd=2)
    legend(2.3, 0.55, c(expression(alpha==0.10), 
                        expression(alpha==0.30),
                        expression(alpha==0.50),
                        expression(alpha==0.75),
                        expression(alpha==1.00),
                        expression(alpha==1.50)), 
           col=1:6, 
           lwd=c(1,1,2,3,3,2), 
           lty=c(1,2,3,4,5,6))

<img src="Apresentacao_Stat4good2_files/figure-markdown_strict/unnamed-chunk-4-1.png" alt="Função de Distribuição Acumulada"  />
<p class="caption">
Função de Distribuição Acumulada
</p>

    #
    t=seq(0,3,by=0.01)
    plot(t,bs(t,0.1,1),main=expression(T%~%BS(alpha,beta==1)),ylab='Densidade',type='l')
    lines(t,bs(t,0.3,1),col=2,lty=2, lwd=1)
    lines(t,bs(t,0.5,1),col=3,lty=3, lwd=2)
    lines(t,bs(t,0.75,1),col=4,lty=4, lwd=3)
    lines(t,bs(t,1,1),col=5,lty=5, lwd=3)
    lines(t,bs(t,1.5,1),col=6,lty=6, lwd=2)
    legend(2.3, 3, c(expression(alpha==0.10), 
                        expression(alpha==0.30),
                        expression(alpha==0.50),
                        expression(alpha==0.75),
                        expression(alpha==1.00),
                        expression(alpha==1.50)), 
           col=1:6, 
           lwd=c(1,1,2,3,3,2), 
           lty=c(1,2,3,4,5,6))

<img src="Apresentacao_Stat4good2_files/figure-markdown_strict/unnamed-chunk-5-1.png" alt="Função Densidade"  />
<p class="caption">
Função Densidade
</p>

<p style="text-align: justify;">
Enquanto que ao manter alpha fixo e mudar o beta há uma mudança na média
e na variança da variável aleatória.
</p>
    t=seq(0.2,3.5,by=0.01)
    plot(t,acbs(t,0.1,0.5),main=expression(T%~%BS(alpha==0.1,beta)),ylab=expression(P(T<=t)==F(t)),type='l')
    lines(t,acbs(t,0.1,0.75),col=2,lty=2, lwd=1)
    lines(t,acbs(t,0.1,1),col=3,lty=3, lwd=2)
    lines(t,acbs(t,0.1,1.5),col=4,lty=4, lwd=3)
    lines(t,acbs(t,0.1,2),col=5,lty=5, lwd=3)
    lines(t,acbs(t,0.1,2.5),col=6,lty=6, lwd=2)
    legend(2.7, 0.55, c(expression(beta==0.50), 
                        expression(beta==0.75),
                        expression(beta==1),
                        expression(beta==1.5),
                        expression(beta==2.00),
                        expression(beta==2.50)), 
           col=1:6, 
           lwd=c(1,1,2,3,3,2), 
           lty=c(1,2,3,4,5,6))

<img src="Apresentacao_Stat4good2_files/figure-markdown_strict/unnamed-chunk-6-1.png" alt="Distribuição Acumulada"  />
<p class="caption">
Distribuição Acumulada
</p>

    t=seq(0,3,by=0.01)
    plot(t,bs(t,0.1,0.5),main=expression(T%~%BS(alpha==0.1,beta)),ylab='Densidade',type='l')
    lines(t,bs(t,0.1,0.75),col=2,lty=2, lwd=1)
    lines(t,bs(t,0.1,1),col=3,lty=3, lwd=2)
    lines(t,bs(t,0.1,1.5),col=4,lty=4, lwd=3)
    lines(t,bs(t,0.1,2),col=5,lty=5, lwd=3)
    lines(t,bs(t,0.1,2.5),col=6,lty=6, lwd=2)
    legend(2, 8, c(expression(beta==0.50), 
                        expression(beta==0.75),
                        expression(beta==1),
                        expression(beta==1.5),
                        expression(beta==2),
                        expression(beta==2.50)), 
           col=1:6, 
           lwd=c(1,1,2,3,3,2), 
           lty=c(1,2,3,4,5,6))

<img src="Apresentacao_Stat4good2_files/figure-markdown_strict/unnamed-chunk-7-1.png" alt="Função Densidade"  />
<p class="caption">
Função Densidade
</p>

Observações:
------------

-   *F*<sub>*T*</sub>(*β*)=0.5, ou seja, *β* é a mediana da
    *B**S*(*α*, *β*)
-   Se *T* ∼ *B**S*(*α*, *β*), então
    *a* &gt; 0,  *a**T* ∼ *B**S*(*α*, *a**β*) e *T*<sup>−1</sup> ∼ *B**S*(*α*, *β*<sup>−1</sup>)
-   A partir da FDA, fazendo
    $Z=\\dfrac{1}{\\alpha}\\left(\\sqrt{\\frac{t}{\\beta}}-\\sqrt{\\frac{\\beta}{t}}\\right)$
    temos:
    \begin{equation}\label{fun3}
    T=\dfrac{\beta}{4}\left[\alpha Z +\sqrt{(\alpha Z)^2+4}\right]^{2}
    \end{equation}
    Essa relação é extremamente útil e pode ser usada para obtenção de
    números pseudo-aleatórios.

a média, a variância, o coeficiente de variação e os coeficientes de
assimetria (*μ*<sub>3</sub>) e curtose (*μ*<sub>4</sub>) da distribuição
BS são, respectivamente:

$$E(T)=\\beta\\left(1+\\dfrac{\\alpha^{2}}{2}\\right),\\ Var(T)=(\\alpha\\beta)^{2}\\left(1+\\dfrac{5\\alpha^{2}}{4}\\right)$$

$$CV(T)=\\dfrac{\\sqrt{5\\alpha^{4}+4\\alpha^{2}}}{\\alpha^{2}+2}$$

$$\\mu\_{3}=\\dfrac{16\\alpha^{2}(11\\alpha^{2}+6)}{(5\\alpha^{2}+4)^{3}},\\ \\mu\_{4}=3+\\dfrac{6\\alpha^{2}(93\\alpha^{2}+41)}{5+\\alpha^{2}+4}$$

A Função Log de Verossimilhança
===============================

Seja *T*<sub>1</sub>, …, *T*<sub>*n*</sub> uma amostra aleatória de
tamanho *n* da distribuição *B**S*(*α*,*β*). Então, o logaritmo da
função de verossimilhança para *θ* = (*α*, *β*), possui a seguinte
forma:
$$
l(\\theta)=-\\dfrac{3}{2}\\sum\_{i=1}^{n}\\log{(t\_{i})}-n\\log{(2\\alpha)}-\\dfrac{n}{2}\\log{(2\\pi)\\beta}-\\dfrac{1}{2\\alpha^{2}}\\sum\_{i=1}^{n}\\left(\\dfrac{t\_{i}}{\\beta}+\\dfrac{\\beta}{t\_{i}}-2\\right)+\\sum\_{i=1}^{n}(t\_{i}+\\beta)$$

Estimadores de Máxima Verossimilhança
-------------------------------------

Os estimadores de máxima verossimilhança (EMV) de *α* e *β* são obtidos
maximizando *l*(*θ*) a partir das soluções das equações:
$$
\\dfrac{\\partial l(\\alpha,\\beta)}{\\partial \\alpha}=-\\dfrac{n}{\\alpha}\\left(1+\\dfrac{2}{\\alpha^{2}}\\right)+\\dfrac{1}{\\beta\\alpha^{3}}\\sum\_{i=1}^{n}t\_{i}+\\dfrac{\\beta}{\\alpha^{3}}\\sum\_{i=1}^{n}\\dfrac{1}{t\_{i}}=0
$$
$$
\\dfrac{\\partial l(\\alpha,\\beta)}{\\partial \\beta}=-\\dfrac{n}{2\\beta}+\\frac{\\displaystyle{\\sum\_{i=1}^{n}t\_{i}}}{2\\alpha^{2}\\beta^{2}}+\\sum\_{i=1}^{n}\\dfrac{1}{t\_{i}+\\beta}-\\dfrac{1}{2\\alpha^{2}}\\sum\_{i=1}^{n}\\dfrac{1}{t\_{i}}=0
$$

<p style="text-align: justify;">
Antes de apresentar a simulação de Monte Carlo para estimação dos
parâmetros veja o gráfico da função log de verossimilhança:
</p>
    ########################################
    require(rgl) # Para fazer os gráficos com rotação
    require(colorRamps) # Para usar a palleta de cor do matlab
    require(plot3D)
    require(rglwidget)

    alpha=2
    beta=1
    n=1000
    set.seed(355)
    z1<-rnorm(1000,0,1)
    dados<-beta*((alpha*z1/2)+sqrt((alpha*z1/2)^2+1))^2
    alpha= seq(1.5, 2.5,l=100)
    beta=seq(0.5,1.5,l=100)
    # A função de verossimilhança para fazer o gráfico:
    f<-function(alpha,beta){
      t=dados
      sum(log(t+beta))-n*log(2*alpha)-(n/2)*(log(2*pi*beta))-(3/2)*sum(log(t))-
        ((1/(2*alpha^2))*sum((t/beta)+(beta/t)-2))
    }

    f <- Vectorize(f)
    z <- outer(alpha,beta, f)

    persp3d(f,xlim=alpha, ylim=beta, contour=TRUE, col = matlab.like,
            xlab = expression(alpha), ylab = expression(beta), zlab = expression(l(Theta)))
    rglwidget()

<img src="Apresentacao_Stat4good2_files/figure-markdown_strict/unnamed-chunk-9-1.png" alt="Função de Velhossimilhança"  />
<p class="caption">
Função de Velhossimilhança
</p>

    hist3D(x=alpha, y=beta, z=z, contour=TRUE, facets=TRUE, curtain=F, phi=-30,theta=30,xlab = expression(alpha), ylab = expression(beta), zlab = expression(l(Theta)))

<img src="Apresentacao_Stat4good2_files/figure-markdown_strict/unnamed-chunk-10-1.png" alt="Histograma da Função de Velhossimilhança"  />
<p class="caption">
Histograma da Função de Velhossimilhança
</p>

    image(x=alpha, y=beta,z,col = terrain.colors(500))
    contour(x=alpha, y=beta, z=z,add=TRUE,levels = pretty(c(-1780,-2000),20),xlab = expression(alpha), ylab = expression(beta), zlab = expression(l(Theta)))
    points(x=c(2,1),pch=19)

<img src="Apresentacao_Stat4good2_files/figure-markdown_strict/unnamed-chunk-11-1.png" alt="Grágico de Contorno da Função de Velhossimilhança"  />
<p class="caption">
Grágico de Contorno da Função de Velhossimilhança
</p>

Possíveis problemas na estimação dos parâmetros
-----------------------------------------------

<p style="text-align: justify;">
Alguns problemas que tive na implementação da verossimilhança e que
merecem destaque:
</p>
<p style="text-align: justify;">
Dois problemas possíveis são aparecer números muito grandes, maiores do
que a maquina possa suportar, ou valores negativos para funções com
suporte positivo, durante a estimação de forma que o R não consiga
continuar o processo de otimização. Vejam o menor número positivo que
pode ser representado pela máquina, o maior número e outros valores
importantes:
</p>
    .Machine

    ## $double.eps
    ## [1] 2,220446e-16
    ## 
    ## $double.neg.eps
    ## [1] 1,110223e-16
    ## 
    ## $double.xmin
    ## [1] 2,225074e-308
    ## 
    ## $double.xmax
    ## [1] 1,797693e+308
    ## 
    ## $double.base
    ## [1] 2
    ## 
    ## $double.digits
    ## [1] 53
    ## 
    ## $double.rounding
    ## [1] 5
    ## 
    ## $double.guard
    ## [1] 0
    ## 
    ## $double.ulp.digits
    ## [1] -52
    ## 
    ## $double.neg.ulp.digits
    ## [1] -53
    ## 
    ## $double.exponent
    ## [1] 11
    ## 
    ## $double.min.exp
    ## [1] -1022
    ## 
    ## $double.max.exp
    ## [1] 1024
    ## 
    ## $integer.max
    ## [1] 2147483647
    ## 
    ## $sizeof.long
    ## [1] 4
    ## 
    ## $sizeof.longlong
    ## [1] 8
    ## 
    ## $sizeof.longdouble
    ## [1] 16
    ## 
    ## $sizeof.pointer
    ## [1] 8

<p style="text-align: justify;">
Exemplo de um número menor que o epsilon da maquina.
</p>
    0.3-0.1==0.2

    ## [1] FALSE

    isTRUE(0.3-0.1==0.2)

    ## [1] FALSE

    0.2-(0.3-0.1)

    ## [1] 2,775558e-17

    print(0.3-0.1,digits=17)

    ## [1] 0,19999999999999998

<p style="text-align: justify;">
O método de otimização que uso (função 'optim' e 'nlimb' do 'R') são
métodos irrestritos, como vários outros métodos de otimização
implementados no R, como os parâmetros da Birnbaum-Saunders possuem
suporte positivo uma modificação (reparametrização) pode ser útil. No
meu script ao definir a função log de verossimilhança mudei alpha e beta
para exp(lalpha) e exp(lbeta), respectivamente, onde lalpha=log(alpha) e
lbeta=log(beta). Ao fazer esta mudança não altera-se em nada a função
pois a substituição usa exp(log(alpha))=alpha e exp(log(beta))=beta.
</p>
<p style="text-align: justify;">
Como mencionado essa modificação é um artifício para o processo de
otimização tornar-se irrestrito, porém deve-se tomar certos cuidados,
uma vez que, valores de exp(*x*)&gt;1.797693*e* + 308 não podem ser
representados pelos computadores usuais. Note que, *x* = 709.7827
(log(1.797693*e* + 308)) é o maior valor de *x* que pode ser
representado nesses computadores, que é relativamente baixo, dado os
problemas que observa-se em inferência estatística. Notem a importância
da função logaritmo em processos de otimização!
</p>
<p style="text-align: justify;">
Caso o método de otimização funcione, mas apareça alguns 'Warnings' com
valores NA, podemos manter alpha e beta e acrescentar na função optim o
comando abaixo:
</p>
`1-suppressWarnings(optim(start,fn=Loglik,method="BFGS",hessian=T)$par)`
<p style="text-align: justify;">
Neste caso escondemos os avisos. Os mesmos continuarão lá, porém
ocultos, não é recomendado! Caso na simulação apareça números negativos
no argumento do log, podemos acrescentar também o código:
</p>
`if (any(c(t, beta, alpha, pi) < 0)) return(NA)`
<p style="text-align: justify;">
Vejamos como gerar valores aleatório com distribuição Birnbaum-Saunders
e a estimação de alpha e beta usando uma simulação de Monte Carlo e as
funções optim e nlimb.
</p>
Simulação de Monte Carlo
========================

<p style="text-align: justify;">
O método de Monte Carlo é um método de simulação estatística que utiliza
sequencias de números aleatórios para desenvolver simulações. Em outras
palavras, é visto como método numérico universal para resolver problemas
por meio de amostragem aleatória.
</p>
    #geração de valores de uma distribuição Birnbaum-Saunders(alpha,beta)

    n<-100

    #Gerando uma variável aleatória t com distribuição Birnbaum-Saunders(alpha,beta)
    alpha=1
    beta=2
    truevalue=c(alpha,beta)
    #Note que teremos N amostras de tamanho 100 distintas, portanto t
    #deve estar dentro do loop do Monte Carlo, alpha e beta deve estar fora do
    #Loop pois os mesmos são fixos para todas as amostras.
    N=1000
     m=matrix(nrow=N,ncol=2)
    m1=matrix(nrow=N,ncol=2)
    for(i in 1:N){
    z<-rnorm(n,0,1)
    #Gerando uma variável aleatória t com distribuição Birnbaum-Saunders(alpha,beta)
    t<-cbind(beta*((alpha*z/2)+sqrt((alpha*z/2)^2+1))^2)
    #Verossimilhança
    Loglik<-function(par,dados){
        lalpha=par[1]
        lbeta=par[2]
        t<-dados
        ll<-sum(log(t+exp(lbeta)))-n*log(2*exp(lalpha))-(n/2)*(log(2*pi*exp(lbeta)))-(3/2)*
          sum(log(t))-((1/(2*exp(lalpha)^2))*sum((t/exp(lbeta))+(exp(lbeta)/t)-2))
        return(-ll)
    }

    #Utilizando a verossimilhança e a função optim para estimar os parâmetros alpha e beta que
    #deram origem as observações t observadas. Note que foi utilizado o log de alpha e 
    #beta como chutes iniciais.
      
    lalpha_0=log(2)
    lbeta_0=log(2)
    start=c(lalpha_0,lbeta_0)
      
    m[i,]=exp(optim(start,fn=Loglik,method="BFGS",dados=t,hessian=T)$par)
    m1[i,]=exp(nlminb(start,Loglik,dados=t)$par)
      

    }

    #Calculating the average of each column of the array of parameters m
    mest=colMeans(m)
    mest1=colMeans(m1)

    #calculating the standard deviation of each column of the array of parameters m
    dest=apply(m,2,sd)
    dest1=apply(m1,2,sd)

    #root mean square error in the calculation of each column of the array of parameters m in 
    #relation to the true value of the parameter
    eqm=function(x,poisson_opt){ 
    k=length(x)
    sqrt(sum(((x-poisson_opt)^2))/k)}

    eqm1=function(x,poisson_nlm){ 
    k=length(x)
    sqrt(sum(((x-poisson_nlm)^2))/k)}

    #Estimated mean squared error of each parameter 
    eqmest=c(eqm(x=m[,1],poisson_opt=truevalue[1]),
             eqm(x=m[,2],poisson_opt=truevalue[2]))

    #Estimated mean squared error of each parameter
    eqmest1=c(eqm1(x=m1[,1],poisson_nlm=truevalue[1]),
              eqm1(x=m1[,2],poisson_nlm=truevalue[2]))

    # Table with the true values of the parameters and the average
    # Standard deviation and mean square error of the estimated parameters
    tab=data.frame(truevalue,mean=mest,sd=dest,eqm=eqmest)
    tab1=data.frame(truevalue,mean=mest1,sd=dest1,eqm=eqmest1)

    tab

    ##   truevalue      mean         sd        eqm
    ## 1         1 0,9890158 0,06917214 0,07000467
    ## 2         2 2,0083009 0,18420781 0,18430271

    tab1

    ##   truevalue      mean        sd       eqm
    ## 1         1 0,9890151 0,0691906 0,0700230
    ## 2         2 2,0083120 0,1842361 0,1843315

    par(mfrow=c(1,2))
    hist(m[,1],prob=T,ylab="densidade",main="");
    rug(m[,1])
    curve(expr = dnorm(x,mean=mean(m[,1]),sd=sd(m[,1])),add=T, col="red")
    hist(m[,2],prob=T,ylab="densidade",main="");
    rug(m[,2])
    curve(expr = dnorm(x,mean=mean(m[,2]),sd=sd(m[,2])),add=T, col="red")

<img src="Apresentacao_Stat4good2_files/figure-markdown_strict/unnamed-chunk-15-1.png" alt="Histograma das estimativas dos parâmetros usando a função optim"  />
<p class="caption">
Histograma das estimativas dos parâmetros usando a função optim
</p>

    par(mfrow=c(1,2))
    hist(m1[,1],prob=T,ylab="densidade",main="");
    rug(m1[,1])
    curve(expr = dnorm(x,mean=mean(m1[,1]),sd=sd(m1[,1])),add=T, col="red")
    hist(m1[,2],prob=T,ylab="densidade",main="");
    rug(m1[,2])
    curve(expr = dnorm(x,mean=mean(m1[,2]),sd=sd(m1[,2])),add=T, col="red")

<img src="Apresentacao_Stat4good2_files/figure-markdown_strict/unnamed-chunk-16-1.png" alt="Histograma das estimativas dos parâmetros usando a função nlimb"  />
<p class="caption">
Histograma das estimativas dos parâmetros usando a função nlimb
</p>

Uma Reparametrização Importante
===============================

<p style="text-align: justify;">
As vezes, reparametrizações são essenciais, pois facilitam o
desenvolvimento analítico de algumas distribuições e também podem
melhorar a eficiência em simulações, em determinadas situações, como em
regressão, quando a distribuição da variável resposta não possui a média
como um de seus parâmetros podemos proceder uma reparametrização de
forma a atender essa condição e poder ajustar a média da variável
resposta. Exemplos de distribuições em que são realizadas
reparametrizações com sucesso são: distribuição beta (ver Ferrari e
Cribari-Neto 2004) e distribuição gaussiana inversa (ver Tweedie 1957).
</p>
Seja $\\mu=\\beta(1+\\dfrac{\\alpha^{2}}{2})$ e
$\\phi=\\dfrac{2}{\\alpha^{2}}.$ Então,
\begin{equation}\label{3}
\alpha=\sqrt{\dfrac{2}{\phi}}\ \textrm{e}\ \beta=\dfrac{\mu}{(1+\frac{1}{\phi})}.
\end{equation}
A fda da *B**S*(*μ*, *ϕ*) é obtida substituindo os valores de *α* e *β*
definidos em $(\\ref{3})$ na expressão $(\\ref{fun1}).$ De onde, temos:
\begin{equation}
F(t;\mu,\phi)=P(T\leq t)=\Phi\left[\sqrt{\dfrac{\phi}{2}}\left(\sqrt{\dfrac{(\phi+1)t}{\phi\mu}}-\sqrt{\dfrac{\phi\mu}{(\phi+1)t}}\right)\right],
\end{equation}
onde *ϕ* &gt; 0,  *μ* &gt; 0,  *t* &gt; 0.

Segue que a fdp é dada por:

$$f(t;\\phi,\\mu)=\\dfrac{exp(\\frac{\\phi}{2})\\sqrt{\\phi+1}}{4\\sqrt{\\pi\\mu}}t^{-\\frac{3}{2}}\\left\[t+\\dfrac{\\phi\\mu}{\\phi+1}\\right\]exp\\left\\{-\\dfrac{\\phi}{4}\\left(\\dfrac{t(\\phi+1)}{\\phi\\mu}+\\dfrac{\\phi\\mu}{t(\\phi+1)}\\right)\\right\\}$$

A nova média e variância são:
\begin{equation}
E(T)=\mu,\ \textrm{e}\ Var(T)=\dfrac{g(\mu)}{h(\phi)} 
\end{equation}
A distribuição *B**S*(*μ*, *ϕ*) satisfaz a propriedade de escala e
também satisfaz a propriedade recíproca.

    bs<-function(t,mu,phi){
      fdp=((exp(phi/2)*sqrt(phi+1))/(4*sqrt(pi*mu)*t^(3/2)))*(t+((phi*mu)/(phi+1)))*
        exp((-phi/4)*((t*(phi+1)/(phi*mu))+(phi*mu/(t*(phi+1)))))
      return(fdp)
    }

    #Teste para ver se a integral da densidade é igual a 1.
    mu=1
    phi=2
    integrand <- function(t){((exp(phi/2)*sqrt(phi+1))/(4*sqrt(pi*mu)*t^(3/2)))*
        (t+((phi*mu)/(phi+1)))*exp((-phi/4)*((t*(phi+1)/(phi*mu))+(phi*mu/(t*(phi+1)))))}
    integrate(integrand, lower = 0, upper = Inf)

    ## 1 with absolute error < 7,8e-05

<p style="text-align: justify;">
Note que com a mudança de *μ*, mantendo *ϕ* fixo, há uma alteração na
curtose da variável aleatória e com o aumento de *ϕ* há um aumento da
assimetria e uma diminuição da variância do gráfico. Notem também que,
ao manter *μ* fixo e aumentar *ϕ* a variância de T tende a zero.
</p>
    #
    t=seq(0.5,6,by=0.01)
    plot(t,bs(t,1,100),main=expression(T%~%BS(mu,phi==100)),ylab='Densidade',type='l')
    lines(t,bs(t,1.5,100),col=2,lty=2, lwd=1)
    lines(t,bs(t,2,100),col=3,lty=3, lwd=2)
    lines(t,bs(t,2.5,100),col=4,lty=4, lwd=3)
    lines(t,bs(t,3,100),col=5,lty=5, lwd=3)
    lines(t,bs(t,3.5,100),col=6,lty=6, lwd=2)
    legend(4.3,2.7, c(expression(mu==1), 
                        expression(mu==1.5),
                        expression(mu==2),
                        expression(mu==2.5),
                        expression(mu==3),
                        expression(mu==3.5)), 
           col=1:6, 
           lwd=c(1,1,2,3,3,2), 
           lty=c(1,2,3,4,5,6))

<img src="Apresentacao_Stat4good2_files/figure-markdown_strict/unnamed-chunk-18-1.png" alt="Função Densidade"  />
<p class="caption">
Função Densidade
</p>

    t=seq(0,4,by=0.01)
    plot(t,bs(t,1,2),ylim=c(0,2.8),main=expression(T%~%BS(mu==1,phi)),ylab='Densidade',type='l')
    lines(t,bs(t,1,5),col=2,lty=2, lwd=1)
    lines(t,bs(t,1,10),col=3,lty=3, lwd=2)
    lines(t,bs(t,1,25),col=4,lty=4, lwd=3)
    lines(t,bs(t,1,50),col=5,lty=5, lwd=3)
    lines(t,bs(t,1,100),col=6,lty=6, lwd=2)
    legend(2.5, 2.5, c(expression(phi==2), 
                        expression(phi==5),
                        expression(phi==10),
                        expression(phi==25),
                        expression(phi==50),
                        expression(phi==100)), 
           col=1:6, 
           lwd=c(1,1,2,3,3,2), 
           lty=c(1,2,3,4,5,6))

<img src="Apresentacao_Stat4good2_files/figure-markdown_strict/unnamed-chunk-19-1.png" alt="Função Densidade"  />
<p class="caption">
Função Densidade
</p>

    V=function(mu,phi){
      vt=2*(mu^2)*(2*phi+5)/((phi+1)^2)
      return(vt)
    }

    # 
    phi=seq(0.001,100,by=0.01)
    plot(phi,V(0.6,phi),ylim=c(0,2.8),main=expression(T%~%BS(mu==1,phi)),ylab='Var(T)',type='l')
    lines(phi,V(1,phi),col=2,lty=2, lwd=1)
    lines(phi,V(1.5,phi),col=3,lty=3, lwd=2)
    lines(phi,V(2,phi),col=4,lty=4, lwd=3)
    lines(phi,V(2.5,phi),col=5,lty=5, lwd=3)
    lines(phi,V(3,phi),col=6,lty=6, lwd=2)
    legend(60, 2.5, c(expression(alpha==0.6), 
                        expression(alpha==1),
                        expression(alpha==1.5),
                        expression(alpha==2),
                        expression(alpha==2.5),
                        expression(alpha==3)), 
           col=1:6, 
           lwd=c(1,1,2,3,3,2), 
           lty=c(1,2,3,4,5,6))

<img src="Apresentacao_Stat4good2_files/figure-markdown_strict/unnamed-chunk-20-1.png" alt="Variância de T"  />
<p class="caption">
Variância de T
</p>

Nova Função Log de Verossimilhança
----------------------------------

$$l(\\mu,\\phi;\\textbf{T})=\\dfrac{n\\phi}{2}+\\dfrac{n}{2}\\log{(\\phi+1)}-\\dfrac{3}{2}\\sum\_{i=1}^{n}\\log{t\_{i}}-n\\log{(4\\sqrt{\\pi\\mu})}+\\sum\_{i=1}^{n}\\log{\\left\[t\_{i}+\\dfrac{\\phi\\mu}{\\phi+1}\\right\]}-\\dfrac{\\phi}{4}\\sum\_{i=1}^{n}\\left\[\\dfrac{t\_{i}(\\phi+1)}{\\phi\\mu}+\\dfrac{\\phi\\mu}{t\_{i}(\\phi+1)}\\right\]$$

<p style="text-align: justify;">
Os estimadores (ou estimativas) de máxima verossimilhança de *μ* e *ϕ*
são obtidos maximizando essa função, a partir da solução das equações
formadas com as derivadas parciais em relação *μ* e *ϕ*. É possível
mostrar que não é possível obter uma solução analítica para os
estimadores de máxima verossimilhança e, portanto, métodos iterativos de
otimização são utilizados.
</p>
Modelos de Regressão Birnbaum-Saunders
======================================

Modelo de regressão Birnbaum-Saunders para *μ*
----------------------------------------------

<p style="text-align: justify;">
Criou-se agora uma estrutura de regressão para a média da distribuição
*B**S*(*α*, *β*) fazendo
*g*(*μ*)=*β*<sub>0</sub> + *β*<sub>1</sub> \* *X*
</p>
    rm(list=ls())
    cat("\014")



    N=1000
    #m e m1 são as matrizes que receberão as estimativas no final do processo de estimação
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

    z<-cbind(rnorm(n,0,1))
    t<-((phi*mu)/(phi+1))*((z/sqrt(2*phi))+(sqrt((z/sqrt(2*phi))^2+1)))^2
      
    Loglik<-function(beta,dados){
     p=ncol(X)
     mu=exp(X%*%beta[1:p])
     phi=phi
     lv=sum((phi/2)-(3/2)*log(t)+(1/2)*log(phi+1)-log(4*sqrt(pi*mu))+log(t+(phi*mu/(phi+1)))
            -(phi/4)*((t*(phi+1)/(phi*mu))+(phi*mu/(t*(phi+1)))))
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

    #root mean square error in the calculation of each column of the array of parameters m in 
    #relation to the true value of the parameter
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

    ##   truevalue      mean        sd       eqm
    ## 1         2 1,9989342 0,1871686 0,1870780
    ## 2         1 0,9897028 0,3373453 0,3373338

    tab1

    ##   truevalue      mean        sd       eqm
    ## 1         2 1,9989497 0,1871741 0,1870835
    ## 2         1 0,9896718 0,3373500 0,3373394

    par(mfrow=c(1,2))
    hist(m[,1],prob=T,ylab="densidade",main="");
    rug(m[,1])
    curve(expr = dnorm(x,mean=mean(m[,1]),sd=sd(m[,1])),add=T, col="red")
    hist(m[,2],prob=T,ylab="densidade",main="");
    rug(m[,2])
    curve(expr = dnorm(x,mean=mean(m[,2]),sd=sd(m[,2])),add=T, col="red")

<img src="Apresentacao_Stat4good2_files/figure-markdown_strict/unnamed-chunk-22-1.png" alt="Histograma das estimativas dos parâmetros obtidas com a função optim"  />
<p class="caption">
Histograma das estimativas dos parâmetros obtidas com a função optim
</p>

    par(mfrow=c(1,2))
    hist(m1[,1],prob=T,ylab="densidade",main="");
    rug(m1[,1])
    curve(expr = dnorm(x,mean=mean(m1[,1]),sd=sd(m1[,1])),add=T, col="red")
    hist(m1[,2],prob=T,ylab="densidade",main="");
    rug(m1[,2])
    curve(expr = dnorm(x,mean=mean(m1[,2]),sd=sd(m1[,2])),add=T, col="red")

<img src="Apresentacao_Stat4good2_files/figure-markdown_strict/unnamed-chunk-23-1.png" alt="Histograma das estimativas dos parâmetros obtidas com a função nlimb"  />
<p class="caption">
Histograma das estimativas dos parâmetros obtidas com a função nlimb
</p>

Modelo de regressão Birnbaum-Saunders para *μ* e *ϕ*
----------------------------------------------------

<p style="text-align: justify;">
Agora, vamos criar uma estrutura de regressão para *μ* e *ϕ*, de tal
forma que *g*(*μ*)=*β*<sub>0</sub> + *β*<sub>1</sub>*X* e
*h*(*ϕ*)=*α*<sub>0</sub> + *α*<sub>1</sub>*Z*
</p>
    N=1000
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

    remove(beta,beta0,beta1,alpha,alpha0,alpha1)
    #Monte Carlo para estimação dos parâmetros beta0, beta1 e beta2

    for (i in 1:N){
      
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
        lv=sum((phi/2)-(3/2)*log(t)+(1/2)*log(phi+1)-log(4*sqrt(pi*mu))+log(t+(phi*mu/(phi+1)))
               -(phi/4)*((t*(phi+1)/(phi*mu))+(phi*mu/(t*(phi+1)))))
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

    #root mean square error in the calculation of each column of the array of parameters m
    #in relation to the true value of the parameter
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

    ##   truevalue      mean         sd        eqm
    ## 1         2  2,001007 0,04273662 0,04272713
    ## 2        -1 -1,000395 0,07786891 0,07783096
    ## 3         3  3,021014 0,30496013 0,30553114
    ## 4         1  1,025665 0,48990768 0,49033482

    tab1

    ##   truevalue      mean         sd        eqm
    ## 1         2  2,001007 0,04273613 0,04272663
    ## 2        -1 -1,000396 0,07786908 0,07783114
    ## 3         3  3,020959 0,30493394 0,30550124
    ## 4         1  1,025743 0,48985955 0,49029085

    par(mfrow=c(2,2))
    hist(m[,1],prob=T,ylab="densidade",main="");
    rug(m[,1])
    curve(expr = dnorm(x,mean=mean(m[,1]),sd=sd(m[,1])),add=T, col="red")
    hist(m[,2],prob=T,ylab="densidade",main="");
    rug(m[,2])
    curve(expr = dnorm(x,mean=mean(m[,2]),sd=sd(m[,2])),add=T, col="red")
    hist(m[,3],prob=T,ylab="densidade",main="");
    rug(m[,3])
    curve(expr = dnorm(x,mean=mean(m[,3]),sd=sd(m[,3])),add=T, col="red")
    hist(m[,4],prob=T,ylab="densidade",main="");
    rug(m[,4])
    curve(expr = dnorm(x,mean=mean(m[,4]),sd=sd(m[,4])),add=T, col="red")

<img src="Apresentacao_Stat4good2_files/figure-markdown_strict/unnamed-chunk-25-1.png" alt="Histograma das estimativas dos parâmetros obtidas com a função optim"  />
<p class="caption">
Histograma das estimativas dos parâmetros obtidas com a função optim
</p>

    par(mfrow=c(2,2))
    hist(m1[,1],prob=T,ylab="densidade",main="");
    rug(m1[,1])
    curve(expr = dnorm(x,mean=mean(m1[,1]),sd=sd(m1[,1])),add=T, col="red")
    hist(m1[,2],prob=T,ylab="densidade",main="");
    rug(m1[,2])
    curve(expr = dnorm(x,mean=mean(m1[,2]),sd=sd(m1[,2])),add=T, col="red")
    hist(m1[,3],prob=T,ylab="densidade",main="");
    rug(m1[,3])
    curve(expr = dnorm(x,mean=mean(m1[,3]),sd=sd(m1[,3])),add=T, col="red")
    hist(m1[,4],prob=T,ylab="densidade",main="");
    rug(m1[,4])
    curve(expr = dnorm(x,mean=mean(m1[,4]),sd=sd(m1[,4])),add=T, col="red")

<img src="Apresentacao_Stat4good2_files/figure-markdown_strict/unnamed-chunk-26-1.png" alt="Histograma das estimativas dos parâmetros obtidas com a função optim"  />
<p class="caption">
Histograma das estimativas dos parâmetros obtidas com a função optim
</p>

Referências
===========

1-Z.W. Birnbaum & S.C. Saunders. A new family of life distribuitions.
Journal of Applied Probability, 6 (1969), 319-327.

2-B.S, Luis Enrique, Modelos Birnbaum-Saunders bivariados – Campinas,
SP: \[s.n.\], 2014.

3-S.N, Manoel Ferreira, Estimação e Modelagem com a distribuição
Birnbaum-Saunders: Uma nova reparametrização - Recife, PE, 2010.

4-Barros, M., Paula, G. A., Leiva, V. (2009). An R implementation for
generalized Birnbaum-Saunders distributions. Computational Statistics
and Data Analysis, 53(5), 1511-1528.

############################################################################### 

\*Mensagem Final e Agradecimento

Agradeço ao Rumenick Pereira pela valiosa ajuda na confecção desta
apresentação e na confecção dos scripts do ´R´, sem a qual esta
apresentação não ficaria tão completa! Muito obrigado!

"Há três caminhos para o sucesso:

1- Ensinar o que se sabe - Generosidade Mental

2- Praticar o que se ensina - Coerência Ética

3- Perguntar o que se ignora - Humildade Intelectual"

    *Mario Sergio Cortella
