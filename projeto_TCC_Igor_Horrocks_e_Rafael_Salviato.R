# Para inserir sessões no código: Ctrl + Shift + R
# Computador do Rafa ------------------------------------------------------
setwd("C:/Users/rafae/Documents/GitHub/TCC-TimeSeries")
# Computador do Igor ------------------------------------------------------

# 1) Função geradora de séries temporais ----------------------------------
# É preciso fazer uma função que gere as séries e armazene em formato
# de lista. A função deve ter pelo menos três argumentos:
# - Um dataframe com os parâmetros para a série temporal
# - Tamanho da série
# - Número de repetições

geradora_ARMA <- function(parametros,n,m,order_arma="ARMA"){
# 1) parâmetros: data frame com os parâmetros
# 2) n: tamanho de cada série temporal
# 3) m: número de simulações
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
# Vou apenas indicar as sessões do artigo onde os códigos se encontram.

# "3 - The Data"

library(strucchange):library(magrittr)

data("USIncExp")
ts.plot(USIncExp, gpars = list(col = c("black", "red"))) : legend("topleft", 
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

# "4 - Generalized fluctuation tests"

ocus <- efp(ecm.model, type="OLS-CUSUM", data=USIncExp2)
me <- efp(ecm.model, type="ME", data=USIncExp2, h=0.2)

# "4.2 - Boundaries and plotting"

bound.ocus <- boundary(ocus, alpha=0.05)
plot(ocus)
plot(ocus, boundary=FALSE) # It is also possible to suppress the boundaries and add them afterwards e.g. another color
lines(bound.ocus, col=4)
lines(-bound.ocus, col=4)
plot(me, functional =NULL)

# "4.3 - Significance testing with empirical fluctuation processes"

sctest(ocus)
sctest(ecm.model, type="OLS-CUSUM", data=USIncExp2)

# "5.1 - F Statistics: function FStats"

fs <-  Fstats(ecm.model, from = c(1990,1), to = c(1999,6), data = USIncExp2)

# "5.2 - Boundaries and plotting"

plot(fs) : plot(fs, pval=TRUE) : plot(fs, aveF=TRUE)

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















