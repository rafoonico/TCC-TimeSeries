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

# sugestão de alteração
