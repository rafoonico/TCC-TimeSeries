#           |***********************************|
#           |                                   |
#           | Se voc� est� usando o RStudio,    |
#           | d� um Alt+O para ocultar as se��es|
#           | e um Shift+Alt+O para expandir-las|
#           |                                   |
#           |***********************************|


# Para inserir sess�es no c�digo: Ctrl + Shift + R
# Computador do Rafa ------------------------------------------------------
setwd("C:/Users/rafae/Documents/GitHub/TCC-TimeSeries")
# Computador do Igor ------------------------------------------------------

# 1) Fun��o geradora de s�ries temporais ----------------------------------
# � preciso fazer uma fun��o que gere as s�ries e armazene em formato
# de lista. A fun��o deve ter pelo menos tr�s argumentos:
# - Um dataframe com os par�metros para a s�rie temporal
# - Tamanho da s�rie
# - N�mero de repeti��es

geradora_ARMA <- function(parametros,n,m,order_arma="ARMA"){
# 1) par�metros: data frame com os par�metros
# 2) n: tamanho de cada s�rie temporal
# 3) m: n�mero de simula��es
  library(magrittr)
  
  if(order_arma=="ARMA"){
    simulacoes <- vector(mode="list",length = m) %>% lapply(function(y){
      simulacao <- split(parametros,f=1:nrow(parametros)) %>%
        lapply(function(x){arima.sim(list(order=x[1:3],
                                          ar=x[[4]],ma=x[[5]]),
                                     n)})
        names(simulacao) <- sapply(1:nrow(parametros),
                                   function(x){paste(c("p=",parametros[x,1],
                                                       ", q=",parametros[x,2],
                                                       ", d=",parametros[x,3],
                                                       ", AR(",parametros[x,4],
                                                       ") e MA=",parametros[x,5],")"),
                                                     collapse="")})
        return(simulacao)
    }) 
  }else if(order_arma=="AR"){
    simulacoes <- vector(mode="list",length = m) %>% lapply(function(y){
      simulacao <- split(parametros,f=1:nrow(parametros)) %>%
        lapply(function(x){arima.sim(list(order=x[1:3],
                                          ar=x[[4]]),
                                     n)})
        names(simulacao) <- sapply(1:nrow(parametros),
                                   function(x){paste(c("p=",parametros[x,1],
                                                       ", q=",parametros[x,2],
                                                       ", d=",parametros[x,3],
                                                       " e AR(",parametros[x,4],")"),
                                                     collapse="")})
        return(simulacao)
    })
  }else if(order_arma=="MA"){
    simulacoes <- vector(mode="list",length = m) %>% lapply(function(y){
      simulacao <- split(parametros,f=1:nrow(parametros)) %>%
        lapply(function(x){arima.sim(list(order=x[1:3],
                                          ma=x[[4]]),
                                     n)})
        names(simulacao) <- sapply(1:nrow(parametros),
                                   function(x){paste(c("p=",parametros[x,1],
                                                       ", q=",parametros[x,2],
                                                       ", d=",parametros[x,3],
                                                       " e MA(",parametros[x,4],")"),
                                                     collapse="")})
        return(simulacao)
    })
  }
  
  
  names(simulacoes) <- sapply(1:m,
                              function(x){paste("m=",x,collapse="")})
  
  return(simulacoes)
}

parametros_func <- function(p,d,q,AR=NULL,MA=NULL){
  if(!is.null(AR) && is.null(MA)){
    resultado <- expand.grid(p,d,q,AR) ; names(resultado)=c("p","d","q","AR")  
  }else if(is.null(AR) && !is.null(MA)){
    resultado <- expand.grid(p,d,q,MA) ; names(resultado)=c("p","d","q","MA")  
  }else if(!is.null(AR) && !is.null(MA)){
    resultado <- expand.grid(p,d,q,AR,MA) ; names(resultado)=c("p","d","q","AR","MA")
  }
  return(resultado)
}


# 2) Testando o pacote 'strucchange" --------------------------------------
# Vou apenas indicar as sess�es do artigo onde os c�digos se encontram.

# "3 - The Data"

library(strucchange) ; library(magrittr)

data("USIncExp")
ts.plot(USIncExp, gpars = list(col = c("black", "red"))) ; legend("topleft", 
                                                                  legend=c("income","expenditure"),
                                                                  col = c("black", "red"), lty=1)
USIncExp2 <- window(USIncExp,start=c(1985,12))
coint.res <- residuals(lm(expenditure~income,data=USIncExp2)) %>% 
                       ts(start=c(1985,12),freq=12) %>%
                       lag(k=-1)
USIncExp2 <- cbind(USIncExp2, diff(USIncExp2), coint.res)
USIncExp2 <- window(USIncExp2, start=c(1986,1), end=c(2001,2))
colnames(USIncExp2) <- c("income","expenditure","diff.income","diff.expenditure","coint.res")
ecm.model <- diff.expenditure~coint.res+diff.income
summary(lm(ecm.model, USIncExp2))
par(mfrow=c(2,2));plot(lm(ecm.model, USIncExp2));par(mfrow=c(1,1))

# "4 - Generalized fluctuation tests"

# 1.1: RE-CUSUM [Brown, Durbin, and Evans (1975)]
recus <- efp(ecm.model, type="Rec-CUSUM", data=USIncExp2)

# 1.2: OLS-CUSUM [Ploberger and Kramer(1992)]
ocus <- efp(ecm.model, type="OLS-CUSUM", data=USIncExp2)

# 1.3: RE-MOSUM [Chu, Hornik, and Kuan (1995a)]
remos <- efp(ecm.model, type="Rec-MOSUM", data=USIncExp2,h=0.2)

# 1.4: OLS-MOSUM [Chu, Hornik, and Kuan (1995a)]
omos <- efp(ecm.model, type="OLS-MOSUM", data=USIncExp2,h=0.2)

X11();par(mfrow=c(2,2));plot(recus);plot(ocus);plot(remos);plot(omos);par(mfrow=c(1,1))

# 1.5: fluctuation (estimates-based process) [Ploberger, Kramer and Kontrus (1989)] 
 
fluc <- efp(ecm.model, type="fluctuation", data=USIncExp2)

# 1.6: ME [Chu, Hornik, and Kuan (1995b)]

me <- efp(ecm.model, type="ME", data=USIncExp2, h=0.2)
me2 <- efp(ecm.model, type="ME", data=USIncExp2, h=0.5)
me3 <- efp(ecm.model, type="ME", data=USIncExp2, h=0.8)

X11();layout(matrix(c(1,1,1,2,3,4),ncol=3,byrow=TRUE));plot(fluc);plot(me);plot(me2);plot(me3)

# Curiosidade: um MOSUM com h=1 � um CUSUM?

teste1 <- efp(ecm.model, type="OLS-MOSUM", data=USIncExp2,h=1)

teste2 <- efp(ecm.model, type="Rec-MOSUM", data=USIncExp2,h=1)

teste3<- efp(ecm.model, type="ME", data=USIncExp2,h=1)

layout(matrix(1:6,ncol=3,byrow=FALSE));plot(ocus);plot(teste1);plot(recus);plot(teste2);plot(fluc);plot(teste3)

#....�, nao da n�o. E se h=0?

teste4 <- efp(ecm.model, type="OLS-MOSUM", data=USIncExp2,h=0)

teste5 <- efp(ecm.model, type="Rec-MOSUM", data=USIncExp2,h=0)

teste6 <- efp(ecm.model, type="ME", data=USIncExp2,h=0)

layout(matrix(1:6,ncol=3,byrow=FALSE));plot(ocus);plot(teste4);plot(recus);plot(teste5);plot(fluc);plot(teste6)

#....Vish, esquece... nao d� n�o. E se h=0.5?

teste7 <- efp(ecm.model, type="OLS-MOSUM", data=USIncExp2,h=0.5)

teste8 <- efp(ecm.model, type="Rec-MOSUM", data=USIncExp2,h=0.5)

teste9 <- efp(ecm.model, type="ME", data=USIncExp2,h=0.5)

layout(matrix(1:6,ncol=3,byrow=FALSE));plot(ocus);plot(teste7);plot(recus);plot(teste8);plot(fluc);plot(teste9)

# beleza, nada a ver uma coisa com a outra.

# "4.2 - Boundaries and plotting"

# Descobri que as linhas do plot(efp) s�o por default o n�vel de signific�nca a 5%

plot(ocus)
# Voltando ao plot(ocus): temos evid�ncia para uma mudan�a estrutural a partir dos anos 90 e ela volta ao normal aos anos 2000

# E aqui, tem a demonstra��o para mudar as linhas (e.g: de 0.05 para 0.1) e inclusive a cor delas
bound.ocus <- boundary(ocus, alpha=0.05)
plot(ocus, boundary=FALSE) # It is also possible to suppress the boundaries and add them afterwards e.g. another color
lines(bound.ocus, col=4)
lines(-bound.ocus, col=4)

# Agora, sobre os efp's baseados em coeficientes: podemos plotar processos por vari�vel
plot(me, functional =NULL)
plot(me) # T� errado usar assim... a explica��o t� mega confusa mas enfim. Evitar.

# "4.3 - Significance testing with empirical fluctuation processes"

# Talvez aqui esteja o principal instrumento de trabalho
sctest(ocus)
sctest(ecm.model, type="OLS-CUSUM", data=USIncExp2)

# Testar com os outros processos
sctest(recus)
sctest(ocus)
sctest(remos)
sctest(omos)
sctest(fluc)
sctest(me)
# Lembrando: H0 � de n�o ter quebra estrutural
# Para obter o p-valor, basta adicionar um $p.value no final
sctest(me)$p.value

# "5.1 - F Statistics: function FStats" [Chow(1960)]

# O teste F foi feito para detectar apenas UMA quebra estrutural.
# O processo � este: separar a amostra em duas, no ponto aonde
# poss�velmente temos a quebra estrutural. E ai, ajustar um modelo para
# cada particao. O teste F ent�o, vai verificar se os coeficientes s�o
# iguais para os dois modelos. Para isso, ele testa a diferen�a entre os
# res�duos do modelo completo versus os res�duos do modelo repartido.

# Ent�o... o que que o pacote prop�e? Fazer um monte de quebra,
# E gerar um monte de estat�stica F para verificar as quebras de tend�ncia. Uma gambiarra.
# No artigo, ele sugere fazer isso numa regi�o suspeita onde acreditamos ter quebra
# de tend�ncia. Mas tamb�m, n�o tem problema colocar o per�odo todo.
fs <-  Fstats(ecm.model, from = c(1990,1), to = c(1999,6), data = USIncExp2)

# "5.2 - Boundaries and plotting"

# Mesmo esquema que os testes de flutua��o generalizada.
# S� tem uma borda, porque a estat�stica F s�
# pode ser positiva.

plot(fs) # Plotando as estat�sticas F
plot(fs, pval=TRUE) # Plotando os p-valores 
plot(fs, aveF=TRUE) # aqui eu n�o entendi... � sobre a m�dia das estisticas F

# "5.3 - Significance testing with F statistics"

sctest(fs, type="expF")
sctest(ecm.model, type = "expF", from=49, to=162, data=USIncExp2)

# "6 - Monitoring with the generalized fluctuation test"

USIncExp3 <- window(USIncExp2, start = c(1986,1), end =c(1989,12))
me.mefp <- mefp(ecm.model, type = "ME", data = USIncExp3, alpha=0.05)

USIncExp3 <- window(USIncExp2, start = c(1986,1), end =c(1990,12))
me.mefp <- monitor(me.mefp)

USIncExp3 <- window(USIncExp2, start = c(1986,1))
me.mefp <- monitor(me.mefp)

me.mefp

USIncExp3 <- window(USIncExp2, start=c(1986,1), end=c(1989,12))
me.efp <- efp(ecm.model, type="ME", data=USIncExp3, h=0.5)
me.mefp <- mefp(me.efp, alpha=0.05)

USIncExp3 <- window(USIncExp2, start=c(1986,1))
me.mefp <- monitor(me.mefp)
















# 3) Montando o esqueleto do trabalho -------------------------------------
# O objetivo dessa sess�o � montar um primeiro
# esqueleto para o trabalho, validar e se estiver
# tudo ok, montar o TCC

library(strucchange)

# 1 Passo: uma fun��o que gera dados. Tem que gerar duas
# vari�veis uma x que � uma s�rie temporal com tend�ncia, e um y
# que � uma var normal onde a m�dia � y=beta0+beta1x.

# 2 Passo: uma fun��o que calcula o "sctest()" para calcular todos os 6 testes.
# puxar o p-valor usando o "sctest()$p.value".

# 3 Passo: uma fun��o que verifica se o p-valor � maior que 0.05 ou n�o,
# e colocar um 0 ou 1. Ai, fazer uma tabelinha com as propor��es geradas.

