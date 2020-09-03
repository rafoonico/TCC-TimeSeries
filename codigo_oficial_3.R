#           |***********************************|
#           |                                   |
#           | Se você está usando o RStudio,    |
#           | dê um Alt+O para ocultar as seções|
#           | e um Shift+Alt+O para expandir-las|
#           |                                   |
#           |***********************************|


# Para inserir sessões no código: Ctrl + Shift + R
# Caminho -----------------------------------------------------------------
# setwd("C:/Users/rafae/Documents/GitHub/TCC-TimeSeries")
# setwd("C:/Users/Igor/OneDrive/Documentos/TCC/Códigos Adicionais")

# 0) Exemplo utilizado na introdução --------------------------------------
# Primeira situacao: quebra não tão visível

n = 100
t = seq(1,n)

# Antes da quebra: vetor de coeficientes é [1 0]. 
# Depois da quebra, vetor de coeficientes é [1 1].
beta = matrix(c(rep(1,n), c(rep(0,n/2),rep(1,n/2))), ncol=2)

x = matrix(c(rep(1,n),t/n), ncol=2)
set.seed(123); e = rnorm(n)
y = e
for(i in 1:n) y[i] = x[i,]%*%beta[i,]+e[i]

par(mar=c(3,3,1,1),mgp=c(1.6,.6,0))
plot(t,y, col=c(rep("red",n/2),rep("blue",n/2)), pch=19)


# Segunda situacao: quebra visível

# Antes da quebra: vetor de coeficientes é [1 0]. 
# Depois da quebra, vetor de coeficientes é [1 10].

beta_2=matrix(c(rep(1,n), c(rep(0,n/2),rep(10,n/2))), ncol=2)

y = e
for(i in 1:n) y[i] = x[i,]%*%beta_2[i,]+e[i]

par(mar=c(3,3,1,1),mgp=c(1.6,.6,0))
plot(t,y, col=c(rep("red",n/2),rep("blue",n/2)), pch=19)


# 1) Montando o esqueleto do trabalho (versão oficial)---------------------
# O objetivo dessa sessão é montar um primeiro
# esqueleto para o trabalho, validar e se estiver
# tudo ok, montar o TCC
rm(list=ls())
library(strucchange)
library(magrittr)

# 1 Passo: uma função que gera dados. Tem que gerar duas
# variáveis uma x que é uma sequência de números, e um y
# que é uma var normal onde a média é y=beta0+beta1x.

# Parâmetros:
# m <- número de simulações para cada cenário.
# n <- tamanho da amostra.
# c <- quantil onde terá a quebra de tendência da amostra. Um número entre 0 e 1. 
# Aí multiplica pelo "n", arredonda para baixo e temos a posição da quebra.
# r_2 <-  Coef de determinacao. Assim controlamos o quão forte é a relação entre
# "x" e "y".
# beta0, beta1 e beta2: parâmetros da relação entre y e x.

geradora_amostra <- function(m,n,c,r_quad,beta0,beta1,beta2){
  # Calculando numero de núcleos
  no_cores <- parallel::detectCores() - 1
  
  # Iniciando cluster
  cl <- parallel::makeCluster(no_cores)
  
  x <- 1:n
  x_c <- x[floor(c*n)]
  y_estimado <- beta0+beta1*x+beta2*ifelse(x<=x_c,0,x-x_c)
  desvpad_y <- sqrt(sum((y_estimado-mean(y_estimado))^2)/(r_quad*n))
  
  retorno <- parallel::parLapply(cl,
                                 vector(mode="list",
                                        length = m),
                                 function(z){
                                   y <- rnorm(n,mean=y_estimado,sd=desvpad_y*sqrt(1-r_quad))
                                   return(data.frame(x,y))
                                 })
  
  parallel::stopCluster(cl)
  
  return(retorno)
}


geradora_amostra2 <- function(m,n,c,r_quad,b0,b1,incl){
  
  
  # Cria marcação do ponto de corte
  x = 1:n
  x_c <- x[floor(c*n)]
  
  if(c!=1){
    # criando a nova reta a partir do ângulo da primeira
    alpha1 = rCAT::rad2deg(atan(b1)) ; if(incl*alpha1 >= 90 |incl*alpha1 <= -90 ) stop("Inclinação muito acentuada -> novo angulo maior que 90º")
    alpha2 = rCAT::deg2rad(incl*alpha1) # Novo ângulo com efeito da inclinação desejada
    b1_new = tan(alpha2) # Valor do coeficiente angular da nova reta
    b2 = b1_new - b1 # Ajuste no valor do coeficiente angular para obter y_estimado depois da quebra
  }else{
    b2=0  
  }
  
  #b0_new = b0 - b2*x_c
  
  y_estimado <- b0 + b1*x + b2*ifelse(x<=x_c,0,x-x_c)
  desvpad_y <- sqrt(sum((y_estimado-mean(y_estimado))^2)/(r_quad*n))
  
  # # Calculando numero de núcleos
  # no_cores <- parallel::detectCores() - 1
  # 
  # # Iniciando cluster
  # cl <- parallel::makeCluster(no_cores)
  
  # retorno <- parallel::parLapply(cl,
  #                                vector(mode="list",
  #                                       length = m),
  #                                function(z){
  #                                  y <- rnorm(n,mean=y_estimado,sd=desvpad_y*sqrt(1-r_quad))
  #                                  return(data.frame(x,y))
  #                                })
  
  retorno <- lapply(vector(mode="list",
                                        length = m),
                                 function(z){
                                   y <- rnorm(n,mean=y_estimado,sd=desvpad_y*sqrt(1-r_quad))
                                   return(data.frame(x,y))
                                 })
  
  return(retorno)
}

# Apenas um teste considerando sem quebra (c=1): teste <- geradora_amostra(1000,1,0.9,10,5,4);summary(lm(y~x,teste))

# 2 Passo: uma função que calcula o "sctest()" para calcular todos os testes e
# puxar o p-valor usando o "sctest()$p.value".

testes <- function(valores_x_e_y, saida="pvalor"){
  if(saida=="pvalor"){
    result <- data.frame(sctest(y~x,type="Rec-CUSUM",data=valores_x_e_y)$p.value,
                         sctest(y~x,type="OLS-CUSUM",data=valores_x_e_y)$p.value,
                         sctest(y~x,type="fluctuation",data=valores_x_e_y)$p.value,
                         sctest(y~x,type="Rec-MOSUM",data=valores_x_e_y,h=0.15)$p.value,
                         sctest(y~x,type="OLS-MOSUM",data=valores_x_e_y,h=0.15)$p.value,
                         sctest(y~x,type="ME",data=valores_x_e_y,h=0.15)$p.value,
                         sctest(y~x,type="supF",data=valores_x_e_y)$p.value,
                         sctest(y~x,type="aveF",data=valores_x_e_y)$p.value,
                         sctest(y~x,type="expF",data=valores_x_e_y)$p.value)
    names(result) <- c("Rec-CUSUM","OLS-CUSUM","fluctuation",
                       "Rec-MOSUM","OLS-MOSUM","ME",
                       "supF","aveF","expF")
  }else if(saida=="T"){
    result <- data.frame(sctest(y~x,type="Rec-CUSUM",data=valores_x_e_y)$statistic,
                         sctest(y~x,type="OLS-CUSUM",data=valores_x_e_y)$statistic,
                         sctest(y~x,type="fluctuation",data=valores_x_e_y)$statistic,
                         sctest(y~x,type="Rec-MOSUM",data=valores_x_e_y,h=0.15)$statistic,
                         sctest(y~x,type="OLS-MOSUM",data=valores_x_e_y,h=0.15)$statistic,
                         sctest(y~x,type="ME",data=valores_x_e_y,h=0.15)$statistic,
                         sctest(y~x,type="supF",data=valores_x_e_y)$statistic,
                         sctest(y~x,type="aveF",data=valores_x_e_y)$statistic,
                         sctest(y~x,type="expF",data=valores_x_e_y)$statistic)
    names(result) <- c("Rec-CUSUM","OLS-CUSUM","fluctuation",
                       "Rec-MOSUM","OLS-MOSUM","ME",
                       "supF","aveF","expF")
  }
  return(result)
}

# 3 Passo: uma função que faz a simulação. Pega os parâmetros, 
# faz o grid, gera as amostras, calcula os testes e obtêm
# as taxas de rejeições da hipótese nula. E também, soltem uma amostra 
# de pelo menos dois cenários só para conferir a distribuição dos p-valores.


# Nova função: uma que mude a inclinação da segunda reta
# a primeira função faz um angulo com o eixo X de 
# 79.11446 graus ou 1.380808 radianos (asin(260/sqrt(50^2+260^2)), onde 260
# é E[Y|X]=10+5*50) "..."

radianos <- asin(260/sqrt(50^2+260^2))
rad2deg <- function(rad) {(rad * 180) / (pi)}
angulo <- rad2deg(radianos)
outro_angulo <- 180-90-angulo

# ... passo anterior está errado. Vamos utilizar a seguinte
# relação: a tangente de um angulo que uma função faz em relação
# ao eixo das abscissas é o coeficiente angular da reta (nosso caso é o beta_1)
# e então, a partir das relacoes feitas numa regressão segmentada
# podemos encontrar o valor de beta 2





simuladora <- function(m,n,c,r_quad,beta0=10,beta1=5,beta2=7){
  
  param <- expand.grid("n"=n,
                       "c"=c,
                       "r2"=r_quad)
  amostras <- apply(param, 1, function(x){geradora_amostra(m=m,
                                                           n=x[1],
                                                           c=x[2],
                                                           r_quad=x[3],
                                                           beta0=beta0,
                                                           beta1=beta1,
                                                           beta2=beta2)}) 
  
  testes_1 <- lapply(amostras,
                     function(y){
                       do.call(rbind.data.frame,
                               lapply(y,testes))
                     }
  ) %>% lapply(function(x){
    t(apply(x,2,function(z){mean(ifelse(z<0.01,1,0))}))
  })
  
  testes_5 <- lapply(amostras,
                     function(y){
                       do.call(rbind.data.frame,
                               lapply(y,testes))
                     }
  ) %>% lapply(function(x){
    t(apply(x,2,function(z){mean(ifelse(z<0.05,1,0))}))
  })
  
  testes_10 <- lapply(amostras,
                      function(y){
                        do.call(rbind.data.frame,
                                lapply(y,testes))
                      }
  ) %>% lapply(function(x){
    t(apply(x,2,function(z){mean(ifelse(z<0.10,1,0))}))
  })
  
  # Amostra de p-valores
  
  indices_amostra_dist <- c(data.frame("index"=1:nrow(param),
                                       param)[(param$n==min(param$n) & 
                                                 param$c==min(param$c) & 
                                                 param$r2==min(param$r2)),1],
                            data.frame("index"=1:nrow(param),
                                       param)[(param$n==max(param$n) & 
                                                 param$c==max(param$c) & 
                                                 param$r2==max(param$r2)),1]) #sample(1:nrow(param),2)
  
  amostra_dist1 <- do.call(rbind.data.frame,
                           lapply(amostras[[indices_amostra_dist[1]]],
                                  testes))
  
  amostra_dist2 <- do.call(rbind.data.frame,
                           lapply(amostras[[indices_amostra_dist[2]]],
                                  testes))
  
  return(list("Nível de significância 0,01" = cbind(param,
                                                    do.call(rbind.data.frame,testes_1)),
              "Nível de significância 0,05" = cbind(param,
                                                    do.call(rbind.data.frame,testes_5)),
              "Nível de significância 0,10" = cbind(param,
                                                    do.call(rbind.data.frame,testes_10)),
              "primeira amostra"=list(indices_amostra_dist[[1]],
                                      amostra_dist1),
              "segunda amostra"=list(indices_amostra_dist[[2]],
                                     amostra_dist2)))
}
# 2) Realizando o estudo de simulação ----------------------------------------

# Vamos fazer um teste, para averiguar quanto tempo leva 
# o cenário com maior número de dados a ser gerado


inicio <- Sys.time()
set.seed(123)
amostra <- simuladora(m = 5000,
                      n = 500,
                      c = 0.75,
                      r_quad = 0.85)
fim <- Sys.time()

fim-inicio # Mano... 1.7 horas, temos que melhorar isso
rm("fim","inicio")

# Bom, vamos gerar as médias.

resultado <- simuladora(m = 10000,
                        n = c(50 , 163 , 275 , 388 , 500),
                        c = 1, #c(0.25,0.5,0.75,1),
                        r_quad = c(0.2,0.5,0.7,0.85,0.95))


saveRDS(resultado,"resultado.RData")

# 3) Avaliando os resultados (i.e: visualizações, análises e etc) ---------

rm(list=ls())
resultados <- readRDS("resultado.RData")

library(magrittr)
library(ggplot2)
library(lemon)

# Preparando os dados (i.e: substituindo os traços e transformando os parâmetros em factor)
uniques_1 <- sapply(colnames(resultados[[1]])[1:3], function(x){unique(resultados[[1]][,x])}) %>%
  lapply(function(x){sort(as.numeric(x))})
uniques_5 <- sapply(colnames(resultados[[2]])[1:3], function(x){unique(resultados[[2]][,x])}) %>%
  lapply(function(x){sort(as.numeric(x))})
uniques_10 <- sapply(colnames(resultados[[3]])[1:3], function(x){unique(resultados[[3]][,x])}) %>%
  lapply(function(x){sort(as.numeric(x))})

for(i in 1:3){
  resultados[[1]][,i]=factor(resultados[[1]][,i],levels = uniques_1[[i]])
  resultados[[2]][,i]=factor(resultados[[2]][,i],levels = uniques_5[[i]])
  resultados[[3]][,i]=factor(resultados[[3]][,i],levels = uniques_10[[i]])
}

old_names <- names(resultados[[1]]) # pega só um... é tudo igual

names(resultados[[1]]) <- c(names(resultados[[1]])[1:3],"Rec_CUSUM","OLS_CUSUM",
                            names(resultados[[1]])[6], "Rec_MOSUM", "OLS_MOSUM",names(resultados[[1]])[9:12])
names(resultados[[2]]) <- c(names(resultados[[2]])[1:3],"Rec_CUSUM","OLS_CUSUM",
                            names(resultados[[2]])[6], "Rec_MOSUM", "OLS_MOSUM",names(resultados[[2]])[9:12])
names(resultados[[3]]) <- c(names(resultados[[3]])[1:3],"Rec_CUSUM","OLS_CUSUM",
                            names(resultados[[3]])[6], "Rec_MOSUM", "OLS_MOSUM",names(resultados[[3]])[9:12])

rm(uniques_1,uniques_5,uniques_10)


# 3.1) Tabelas c = 1 ------------------------------------------------------------

# Primeiro vamos visualizá-las

for(i in 1:3){
  names(resultados[[i]]) <- old_names
}

View(resultados[[1]])
View(resultados[[2]])
View(resultados[[3]])

# Depois, gerar as tabelas em latex

# Corrige formato n
n_new1 = as.integer(resultados[[1]]$n)
n_new2 = as.integer(resultados[[2]]$n)
n_new3 = as.integer(resultados[[3]]$n)
resultados[[1]]$n <- n_new1
resultados[[2]]$n <- n_new2
resultados[[3]]$n <- n_new3


# Padroniza com 5 casas decimais
results_format_alpha1 <- apply(resultados[[1]], 1, function(x){
                  new_value = c(x[1:3],format(round(x[4:12], 5), nsmall = 5))
                  return(new_value)
                })

results_format_alpha5 <- apply(resultados[[2]], 1, function(x){
  new_value = c(x[1:3],format(round(x[4:12], 5), nsmall = 5))
  return(new_value)
})

results_format_alpha10 <- apply(resultados[[3]], 1, function(x){
  new_value = c(x[1:3],format(round(x[4:12], 5), nsmall = 5))
  return(new_value)
})

new_results1 = t(results_format_alpha1)
new_results2 = t(results_format_alpha5)
new_results3 = t(results_format_alpha10)


# Ttira a coluna "c" porque é toda igual
library(xtable)
print(xtable(new_results1[,-2]), include.rownames=FALSE)
print(xtable(new_results2[,-2]), include.rownames=FALSE)
print(xtable(new_results3[,-2]), include.rownames=FALSE)


# 3.2) Visualização gráfica final dos resultados para c = 1 --------------------------

library(lattice)
library(reshape)

# Alfa 0,01

melted_resultado <- melt(resultados[[1]][,-2],id.vars= c("r2", "n"))


colors = colorRampPalette(c("blue", "red"))
nlines = nlevels(as.factor(melted_resultado$n))

xyplot( value ~ r2 | variable, 
        data=melted_resultado, as.table=T, group = n,
        auto.key = list(columns = nlines, lwd=4, cex=1),
        par.settings = list(superpose.symbol = list(col = colors(nlines), pch = 19),
                            superpose.line = list(col = colors(nlines))),
        type="l", 
        lwd = 2, layout=c(3,3),
        xlab = "Coef. Determinação (r2)",
        ylab = "Proporção Rejeição H0",
        yscale.components=function(...){
          yc <- yscale.components.default(...)
          yc$left$labels$labels <-
            sprintf("%s%%",yc$left$labels$at*100) ## convert to strings as pct
          return(yc) 
        }
)



for (i in 1:3) {
  for (j in 1:3) {
    trellis.focus("panel", i, j)
    panel.abline(h=0.01)
    trellis.unfocus()
  }
}

plot.save <- recordPlot()

trellis.device(device="png", filename="xyplot_alpha1.png")
print(plot.save)
dev.off()


# Alfa 0,5

melted_resultado <- melt(resultados[[2]][,-2],id.vars= c("r2", "n"))

xyplot( value ~ r2 | variable, 
        data=melted_resultado, as.table=T, group = n,
        auto.key = list(columns = nlines, lwd=4, cex=1),
        par.settings = list(superpose.symbol = list(col = colors(nlines), pch = 19),
                            superpose.line = list(col = colors(nlines))),
        type="l", 
        lwd = 2, layout=c(3,3),
        xlab = "Coef. Determinação (r2)",
        ylab = "Proporção Rejeição H0",
        yscale.components=function(...){
          yc <- yscale.components.default(...)
          yc$left$labels$labels <-
            sprintf("%s%%",yc$left$labels$at*100) ## convert to strings as pct
          return(yc) 
        }
)



for (i in 1:3) {
  for (j in 1:3) {
    trellis.focus("panel", i, j)
    panel.abline(h=0.05)
    trellis.unfocus()
  }
}

plot.save <- recordPlot()

trellis.device(device="png", filename="xyplot_alpha5.png")
print(plot.save)
dev.off()

# Alfa 0,10

melted_resultado <- melt(resultados[[3]][,-2],id.vars= c("r2", "n"))

xyplot( value ~ r2 | variable, 
        data=melted_resultado, as.table=T, group = n,
        auto.key = list(columns = nlines, lwd=4, cex=1),
        par.settings = list(superpose.symbol = list(col = colors(nlines), pch = 19),
                            superpose.line = list(col = colors(nlines))),
        type="l", 
        lwd = 2, layout=c(3,3),
        xlab = "Coef. Determinação (r2)",
        ylab = "Proporção Rejeição H0",
        yscale.components=function(...){
          yc <- yscale.components.default(...)
          yc$left$labels$labels <-
            sprintf("%s%%",yc$left$labels$at*100) ## convert to strings as pct
          return(yc) 
        }
)



for (i in 1:3) {
  for (j in 1:3) {
    trellis.focus("panel", i, j)
    panel.abline(h=0.1)
    trellis.unfocus()
  }
}

plot.save <- recordPlot()

trellis.device(device="png", filename="xyplot_alpha10.png")
print(plot.save)
dev.off()

# 4) Realizando o estudo do Estimated P-Value (EPV) -----------------------

# Para lidar com os cenários aonde a hipótese nula
# não é verdadeira, é preciso pensar em alguma forma de 
# estimar o poder do teste. E para tanto, tem o EPV.

# Primeiro passo: fazer a função "simuladora 2"

simuladora_2 <- function(m,n,c,r_quad,beta0=10,beta1=5,beta2=7){
  
  param <- expand.grid("n"=n,
                       "r2"=r_quad)
  
  #1) Gerar as amostras c=1
  
  amostras_H0 <- apply(param, 1, function(x){geradora_amostra(m=m,
                                                              n=x[1],
                                                              c=1,
                                                              r_quad=x[2],
                                                              beta0=beta0,
                                                              beta1=beta1,
                                                              beta2=beta2)}
  ) %>% lapply(function(y){
    do.call(rbind.data.frame,
            lapply(y,
                   function(w){
                     testes(w,saida="T")
                   }
            )
    )
  })
  
  #2) Gerar as amostras c<1 e calcular as métricas
  
  # ... vamos certificar que o vetor "c" não tem 1 nele
  
  c <- c[c!=1]
  
  # agora, é criar o loop
  
  resultados_brutos <- sapply(c,
                              function(w){
                                amostras <- apply(param, 1, function(x){geradora_amostra(m=m,
                                                                                         n=x[1],
                                                                                         c=w,
                                                                                         r_quad=x[2],
                                                                                         beta0=beta0,
                                                                                         beta1=beta1,
                                                                                         beta2=beta2)}
                                ) %>% lapply(function(y){
                                  do.call(rbind.data.frame,
                                          lapply(y,
                                                 function(w){
                                                   testes(w,saida="T")
                                                 }
                                          )
                                  )
                                })
                                
                                nome_testes <- names(amostras[[1]])
                                
                                A <- vector(mode="list",length = nrow(param))
                                B <- vector(mode="list",length = nrow(param))                               
                                
                                for(i in 1:nrow(param)){
                                  
                                  a <- vector(mode="list",length = ncol(amostras[[1]]))
                                  b <- vector(mode="list",length = ncol(amostras[[1]]))
                                  for(j in 1:ncol(amostras[[1]])){
                                    
                                    a[[j]] <- format(sum((amostras_H0[[i]][,j]>=amostras[[i]][,j])*1)/m,
                                                     20)
                                    
                                    ingrediente_B <- expand.grid(amostras_H0[[i]][,j],
                                                                 amostras[[i]][,j])
                                    
                                    b[[j]] <- format(sum((ingrediente_B[,1]>=ingrediente_B[,2])*1)/(m^2),
                                                     20)
                                    
                                    rm(ingrediente_B)
                                  }
                                  A[[i]] <- do.call(cbind.data.frame,a); names(A[[i]]) <- nome_testes
                                  B[[i]] <- do.call(cbind.data.frame,b); names(B[[i]]) <- nome_testes
                                  
                                  rm(a,b)
                                }
                                
                                # Fazendo um pré-tratamento nos resultados
                                
                                A <- do.call(rbind.data.frame,A)
                                B <- do.call(rbind.data.frame,B)
                                
                                
                                return(list("Parametro C"=w,
                                            "Métrica A"=A,
                                            "Métrica B"=B))
                                
                              })
  
  # 3) Dá a limpeza final e retornar os dados
  
  TabelaA <- do.call(rbind.data.frame,
                     apply(data.frame(seq(1,length(c)*3,3)),1,
                           function(x){
                             return(data.frame("c"=rep(resultados_brutos[[x]],
                                                       times=nrow(param)),
                                               param,
                                               resultados_brutos[[x+1]]))
                           }))
  TabelaB <- do.call(rbind.data.frame,
                     apply(data.frame(seq(1,length(c)*3,3)),1,
                           function(x){
                             return(data.frame("c"=rep(resultados_brutos[[x]],
                                                       times=nrow(param)),
                                               param,
                                               resultados_brutos[[x+2]]))
                           }))
  
  
  
  return(list("EPV Método A"=TabelaA,
              "EPV Método B"=TabelaB))
}



# 5) Realizando o estudo de simulação do EPV ------------------------------

# Vamos fazer um teste, para averiguar quanto tempo leva 
# o cenário com maior número de dados a ser gerado


inicio <- Sys.time()
set.seed(123)
amostra <- simuladora_2(m = 100,
                        n = 50,
                        c = 0.75,
                        r_quad = 0.85)
fim <- Sys.time()

fim-inicio #Time difference of 7.538328 secs
rm("fim","inicio")

# Bom, vamos gerar as médias.

resultado <- simuladora_2(m = 10000,
                          n = c(50 , 163 , 275 , 388 , 500),
                          c = c(0.25,0.5,0.75),
                          r_quad = c(0.2,0.5,0.7,0.85,0.95))


saveRDS(resultado,"resultado_epv.RData")



# 5.1) Tabelas c < 1 ------------------------------------------------------

# Carrega RData com simulação
resultado_epv <- readRDS("resultado_epv.RData")

n_epv_new1 = as.integer(resultado_epv[[1]]$n)
n_epv_new2 = as.integer(resultado_epv[[2]]$n)
resultado_epv[[1]]$n <- n_epv_new1
resultado_epv[[2]]$n <- n_epv_new2


# Padroniza com 5 casas decimais
results_format_epv1 <- apply(resultado_epv[[1]], 1, function(x){
  new_value_epv = c(x[1:3],format(round(as.numeric(as.character(x))[4:12], 5), nsmall = 5))
  return(new_value_epv)
})

results_format_epv2 <- apply(resultado_epv[[2]], 1, function(x){
  new_value_epv = c(x[1:3],format(round(as.numeric(as.character(x))[4:12], 5), nsmall = 5))
  return(new_value_epv)
})

new_results_epv1 = t(results_format_epv1) ; colnames(new_results_epv1) = colnames(resultado_epv[[1]])
new_results_epv2 = t(results_format_epv2) ; colnames(new_results_epv2) = colnames(resultado_epv[[2]])


# Gera tabelas no formato latex - saparadas por ponto de corte
library(xtable)
# c = 0.25
tab_c25_A = cbind(Métrica = "A", new_results_epv1[which(new_results_epv1[,1] == "0.25"),-1]) 
tab_c25_B = cbind(Métrica = "B", new_results_epv2[which(new_results_epv2[,1] == "0.25"),-1])
print(xtable(rbind(tab_c25_A, tab_c25_B)), include.rownames=FALSE)

# c = 0.5
tab_c50_A = cbind(Métrica = "A", new_results_epv1[which(new_results_epv1[,1] == "0.50"),-1]) 
tab_c50_B = cbind(Métrica = "B", new_results_epv2[which(new_results_epv2[,1] == "0.50"),-1])
print(xtable(rbind(tab_c50_A, tab_c50_B)), include.rownames=FALSE)

# c = 0.75
tab_c75_A = cbind(Métrica = "A", new_results_epv1[which(new_results_epv1[,1] == "0.75"),-1]) 
tab_c75_B = cbind(Métrica = "B", new_results_epv2[which(new_results_epv2[,1] == "0.75"),-1])
print(xtable(rbind(tab_c75_A, tab_c75_B)), include.rownames=FALSE)


# Comparando métricas
# fixa n == 50
n50_metA = cbind(Métrica = "A", new_results_epv1[which(new_results_epv1[,2] == " 50"),-2]) 
n50_metB = cbind(Métrica = "B", new_results_epv2[which(new_results_epv2[,2] == " 50"),-2])
tbl_n50 = rbind(n50_metA, n50_metB)
r2 = tbl_n50[,3]
ct = tbl_n50[,2]
tbl_n50_v2 = tbl_n50[order(r2, ct),]
print(xtable(tbl_n50_v2), include.rownames=FALSE)

# fixa n == 163
n163_metA = cbind(Métrica = "A", new_results_epv1[which(new_results_epv1[,2] == "163"),-2])
n163_metB = cbind(Métrica = "B", new_results_epv2[which(new_results_epv2[,2] == "163"),-2])
tbl_n163 = rbind(n163_metA, n163_metB)
r2 = tbl_n163[,3]
ct = tbl_n163[,2]
tbl_n163_v2 = tbl_n163[order(r2, ct),]
print(xtable(tbl_n163_v2), include.rownames=FALSE)

# fixa n == 275
n275_metA = cbind(Métrica = "A", new_results_epv1[which(new_results_epv1[,2] == "275"),-2])
n275_metB = cbind(Métrica = "B", new_results_epv2[which(new_results_epv2[,2] == "275"),-2])
tbl_n275 = rbind(n275_metA, n275_metB)
r2 = tbl_n275[,3]
ct = tbl_n275[,2]
tbl_n275_v2 = tbl_n275[order(r2, ct),]
print(xtable(tbl_n275_v2), include.rownames=FALSE)

# fixa n == 388
n388_metA = cbind(Métrica = "A", new_results_epv1[which(new_results_epv1[,2] == "388"),-2])
n388_metB = cbind(Métrica = "B", new_results_epv2[which(new_results_epv2[,2] == "388"),-2])
tbl_n388 = rbind(n388_metA, n388_metB)
r2 = tbl_n388[,3]
ct = tbl_n388[,2]
tbl_n388_v2 = tbl_n388[order(r2, ct),]
print(xtable(tbl_n388_v2), include.rownames=FALSE)

# fixa n == 500
n500_metA = cbind(Métrica = "A", new_results_epv1[which(new_results_epv1[,2] == "500"),-2])
n500_metB = cbind(Métrica = "B", new_results_epv2[which(new_results_epv2[,2] == "500"),-2])
tbl_n500 = rbind(n500_metA, n500_metB)
r2 = tbl_n500[,3]
ct = tbl_n500[,2]
tbl_n500_v2 = tbl_n500[order(r2, ct),]
print(xtable(tbl_n500_v2), include.rownames=FALSE)


# Tabelas de campeões dos cenários

# Métrica A
tbl_metA = resultado_epv[[1]]
tbl_metA$value = as.numeric(as.character(tbl_metA$value))

summary_metA = apply(tbl_metA, 1, function(x){
  valor_min = min(x[4:12])
  retorno = c(x[1:3], ifelse((x[4:12] == valor_min) == T, "X", ""))
  return(retorno)
})
View(t(summary_metA))
print(xtable(t(summary_metA)),include.rownames=FALSE)

# Métrica B
tbl_metB = resultado_epv[[2]]
# tbl_metB$value = as.numeric(as.character(tbl_metB$value))

summary_metB = apply(tbl_metB, 1, function(x){
  valor_min = min(x[4:12])
  retorno = c(x[1:3], ifelse((x[4:12] == valor_min) == T, "X", ""))
  return(retorno)
})
View(t(summary_metB))
print(xtable(t(summary_metB)),include.rownames=FALSE)

# soma coluna dos testes campeões para verificar diferença de performance entre as métricas
count_metA = apply(tbl_metA, 1, function(x){
  valor_min = min(x[4:12])
  retorno = ifelse((x[4:12] == valor_min) == T, 1, 0)
  return(retorno)
  })
sum_metA <- apply(t(count_metA), 2, sum)

count_metB = apply(tbl_metB, 1, function(x){
  valor_min = min(x[4:12])
  retorno = ifelse((x[4:12] == valor_min) == T, 1, 0)
  return(retorno)
})
sum_metB <- apply(t(count_metB), 2, sum)

sum_final = rbind(c(Métrica = "A", sum_metA), c(Métrica = "B", sum_metB))
print(xtable(sum_final),include.rownames=FALSE)


# 5.2) Visualização Resultados c < 1 ------------------------------------------

# Carrega pacotes para geração dos gráficos

library(lattice)
library(reshape)

corteList = c(0.25, 0.5, 0.75)
metList = c("A", "B")

for (j in 1:2) {
  
  met <- metList[j]
  
  for (i in 1:3) {
    
  corte = corteList[i]
    
  melted_resultado <- melt(resultado_epv[[j]][which(resultado_epv[[j]][1] == corte),-1],id.vars= c("r2", "n"))
  melted_resultado2 <- cbind(melted_resultado,value2 = as.numeric(as.character(melted_resultado$value)))
  
  colors = colorRampPalette(c("blue", "red"))
  nlines = nlevels(as.factor(melted_resultado$n))
  
  plot.save <- xyplot( value2 ~ r2 | variable, 
                data=melted_resultado2, as.table=T, group = n,
                auto.key = list(columns = nlines, lwd=4, cex=1),
                par.settings = list(superpose.symbol = list(col = colors(nlines), pch = 19),
                                    superpose.line = list(col = colors(nlines))),
                type="l", 
                lwd = 2, layout=c(3,3),
                xlab = "Coef. Determinação (r2)",
                ylab = "P-Valor esperado",
                main = paste("Métrica:", met, "\n", "c =", corte)
        )
  
  
  trellis.device(device="png", filename=paste("xyplot_corte",corte*100,"metrica",met,".png", sep = ""))
  print(plot.save)
  dev.off()
  
  }
}




# 6) Versão final do estudo do EPV ----------------------------------------

# Precisamos fazer diferente o afastamento
# da hipótese nula. Vamos mexer na inclinação

simuladora_3 <- function(m,n,c,r_quad,beta0,beta1,incl){
  
  param <- expand.grid("n"=n,
                       "r2"=r_quad,
                       "incl"=incl)
  
  #1) Gerar as amostras c=1
  #geradora_amostra2 <- function(m,n,c,r_quad,b0,b1,incl)
  amostras_H0 <- apply(param, 1, function(x){geradora_amostra2(m=m,
                                                              n=x[1],
                                                              c=1,
                                                              r_quad=x[2],
                                                              b0=beta0,
                                                              b1=beta1,
                                                              incl=x[3])}
  ) %>% lapply(function(y){
    do.call(rbind.data.frame,
            lapply(y,
                   function(w){
                     testes(w,saida="T")
                   }
            )
    )
  })
  
  #2) Gerar as amostras c<1 e calcular as métricas
  
  # ... vamos certificar que o vetor "c" não tem 1 nele
  
  c <- c[c!=1]
  
  # agora, é criar o loop
  
  resultados_brutos <- sapply(c,
                              function(w){
                                amostras <- apply(param, 1, function(x){geradora_amostra2(m=m,
                                                                                         n=x[1],
                                                                                         c=w,
                                                                                         r_quad=x[2],
                                                                                         b0=beta0,
                                                                                         b1=beta1,
                                                                                         incl=x[3])}
                                ) %>% lapply(function(y){
                                  do.call(rbind.data.frame,
                                          lapply(y,
                                                 function(z){
                                                   testes(z,saida="T")
                                                 }
                                          )
                                  )
                                })
                                
                                nome_testes <- names(amostras[[1]])
                                
                                A <- vector(mode="list",length = nrow(param))
                                B <- vector(mode="list",length = nrow(param))                               
                                
                                for(i in 1:nrow(param)){
                                  
                                  a <- vector(mode="list",length = ncol(amostras[[1]]))
                                  b <- vector(mode="list",length = ncol(amostras[[1]]))
                                  for(j in 1:ncol(amostras[[1]])){
                                    
                                    a[[j]] <- format(sum((amostras_H0[[i]][,j]>=amostras[[i]][,j])*1)/m,
                                                     20)
                                    
                                    ingrediente_B <- expand.grid(amostras_H0[[i]][,j],
                                                                 amostras[[i]][,j])
                                    
                                    b[[j]] <- format(sum((ingrediente_B[,1]>=ingrediente_B[,2])*1)/(m^2),
                                                     20)
                                    
                                    rm(ingrediente_B)
                                  }
                                  A[[i]] <- do.call(cbind.data.frame,a); names(A[[i]]) <- nome_testes
                                  B[[i]] <- do.call(cbind.data.frame,b); names(B[[i]]) <- nome_testes
                                  
                                  rm(a,b)
                                }
                                
                                # Fazendo um pré-tratamento nos resultados
                                
                                A <- do.call(rbind.data.frame,A)
                                B <- do.call(rbind.data.frame,B)
                                
                                
                                return(list("Parametro C"=w,
                                            "Métrica A"=A,
                                            "Métrica B"=B))
                                
                              })
  
  # 3) Dá a limpeza final e retornar os dados
  
  TabelaA <- do.call(rbind.data.frame,
                     apply(data.frame(seq(1,length(c)*3,3)),1, # Foi meio confuso... mas o lance de 3 em 3 é ref a lista dos dados brutos (que é param C, metrica A e netrica B). Nao tem a ver com o df "param".
                           function(x){
                             return(data.frame("c"=rep(resultados_brutos[[x]],
                                                       times=nrow(param)),
                                               param,
                                               resultados_brutos[[x+1]]))
                           }))
  TabelaB <- do.call(rbind.data.frame,
                     apply(data.frame(seq(1,length(c)*3,3)),1,
                           function(x){
                             return(data.frame("c"=rep(resultados_brutos[[x]],
                                                       times=nrow(param)),
                                               param,
                                               resultados_brutos[[x+2]]))
                           }))
  
  
  
  return(list("EPV Método A"=TabelaA,
              "EPV Método B"=TabelaB))
}

set.seed(123)
resultado <- simuladora_3(m = 10000,
             n = c(50 , 163 , 275 , 388 , 500),
             c = c(0.25,0.5,0.75),
             r_quad = c(0.2,0.5,0.7,0.85,0.95),
             beta0=10,
             beta1=1,
             incl=c(0.05,0.10,0.15,0.25,0.50))

saveRDS(resultado,"resultado_epv2.RData")


# 7) Versão final2 do estudo do EPV ---------------------------------------

# Agora sim vamos fazer essa versao final...

# Primeiro, gerar amostra sob H0 para usar como argumento na função

gerar_h0 <- function(m,n,r_quad,beta0,beta1){
  param <- expand.grid("n"=n,
                       "r2"=r_quad)
  
  #1) Gerar as amostras c=1
  #m,n,c,r_quad,beta0,beta1,beta2
  amostras_H0 <- apply(param, 1, function(x){geradora_amostra(m=m,
                                                              n=x[1],
                                                              c=1,
                                                              r_quad=x[2],
                                                              beta0=beta0,
                                                              beta1=beta1,
                                                              beta2=7)} # Vai um número qualquer
  ) %>% lapply(function(y){
    do.call(rbind.data.frame,
            lapply(y,
                   function(w){
                     testes(w,saida="T")
                   }
            )
    )
  })
  return(amostras_H0)
}

# Depois, a função propriamente dita

gerar_EPV <- function(m,n,c,r_quad,b0,b1,incl,amostras_H0){
  param <- expand.grid("n"=n,
                       "r2"=r_quad,
                       "incl"=incl)

  c <- c[c!=1]
  
  # agora, é criar o loop
  
  resultados_brutos <- sapply(c,
                              function(w){
                                amostras <- apply(param, 1, function(x){geradora_amostra2(m=m,
                                                                                         n=x[1],
                                                                                         c=w,
                                                                                         r_quad=x[2],
                                                                                         b0=beta0,
                                                                                         b1=beta1,
                                                                                         incl=x[3])}
                                ) %>% lapply(function(y){
                                  do.call(rbind.data.frame,
                                          lapply(y,
                                                 function(z){
                                                   testes(z,saida="T")
                                                 }
                                          )
                                  )
                                })
                                
                                nome_testes <- names(amostras[[1]])
                                
                                A <- vector(mode="list",length = nrow(param))
                                B <- vector(mode="list",length = nrow(param))                               
                                
                                for(i in 1:nrow(param)){

                                  a <- vector(mode="list",length = ncol(amostras[[1]]))
                                  b <- vector(mode="list",length = ncol(amostras[[1]]))
                                  for(j in 1:ncol(amostras[[1]])){
                                    
                                    a[[j]] <- format(sum((amostras_H0[[i]][,j]>=amostras[[i]][,j])*1)/m,
                                                     20)
                                    
                                    ingrediente_B <- expand.grid(amostras_H0[[i]][,j],
                                                                 amostras[[i]][,j])
                                    
                                    b[[j]] <- format(sum((ingrediente_B[,1]>=ingrediente_B[,2])*1)/(m^2),
                                                     20)                                   
                                    rm(ingrediente_B)
                                  }
                                  A[[i]] <- do.call(cbind.data.frame,a); names(A[[i]]) <- nome_testes
                                  B[[i]] <- do.call(cbind.data.frame,b); names(B[[i]]) <- nome_testes
                                  
                                  rm(a,b)
                                }
                                
                                # Fazendo um pré-tratamento nos resultados
                                
                                A <- do.call(rbind.data.frame,A)
                                B <- do.call(rbind.data.frame,B)
                                
                                return(list("Parametro C"=w,
                                            "Métrica A"=A,
                                            "Métrica B"=B))
                                
                              })
  
  # 3) Dá a limpeza final e retornar os dados
  
  TabelaA <- do.call(rbind.data.frame,
                     apply(data.frame(seq(1,length(c)*3,3)),1, # Foi meio confuso... mas o lance de 3 em 3 é ref a lista dos dados brutos (que é param C, metrica A e netrica B). Nao tem a ver com o df "param".
                           function(x){
                             return(data.frame("c"=rep(resultados_brutos[[x]],
                                                       times=nrow(param)),
                                               param,
                                               resultados_brutos[[x+1]]))
                           }))
  TabelaB <- do.call(rbind.data.frame,
                     apply(data.frame(seq(1,length(c)*3,3)),1,
                           function(x){
                             return(data.frame("c"=rep(resultados_brutos[[x]],
                                                       times=nrow(param)),
                                               param,
                                               resultados_brutos[[x+2]]))
                           }))
  return(list("EPV Método A"=TabelaA,
              "EPV Método B"=TabelaB))  
}

# E terceiro: realizar a simulação

m = 10000
n = c(50 , 163 , 275 , 388 , 500)
c = c(0.25,0.5,0.75)
r_quad = c(0.2,0.5,0.7,0.85,0.95)
beta0=10
beta1=1
incl=c(0.05,0.10,0.15,0.25,0.50)

primeiro_for <- expand.grid(n,r_quad, KEEP.OUT.ATTRS = TRUE, stringsAsFactors = TRUE)

segundo_for <- expand.grid(c,incl, KEEP.OUT.ATTRS = TRUE, stringsAsFactors = TRUE)

for(i in 1:nrow(primeiro_for)){ # Desliguei sem querer
  sim_h0 <- gerar_h0(m,primeiro_for[i,1],primeiro_for[i,2],beta0,beta1)
  for(j in 1:nrow(segundo_for)){
    epv <- gerar_EPV(m,
                     primeiro_for[i,1],
                     segundo_for[j,1],
                     primeiro_for[i,2],
                     beta0,
                     beta1,
                     segundo_for[j,2],
                     sim_h0)
    saveRDS(epv,paste(c("resultado_epv",primeiro_for[i,],segundo_for[j,],".RData"),collapse="_"))
    rm(epv)
  }
  rm(sim_h0)
}



# 8) Resultados do EPV ----------------------------------------------------
library(magrittr)
setwd("C:/Users/rafae/Documents/GitHub/TCC-TimeSeries/Bases TCC/versaochega")
resultado <- apply(expand.grid("n" = c(50 , 163 , 275 , 388 , 500),
            "r_quad" = c(0.2,0.5,0.7,0.85,0.95),
            "c" = c(0.25,0.5,0.75),
            "incl"=c(0.05,0.10,0.15,0.25,0.50)),1,function(x){
              paste(c("resultado_epv",x[1],x[2],x[3],x[4],".RData"),collapse="_")
            }) %>% sapply(readRDS)

metodo_A <- t(sapply(seq(1,749,by=2),function(x){resultado[[x]]}))
metodo_B <- t(sapply(seq(2,750,by=2),function(x){resultado[[x]]}))

saveRDS(metodo_A,"metodo_A.RData")
saveRDS(metodo_B,"metodo_B.RData")


# Os resultados estão estranhos...

# 1) Avaliar os betas para diferentes inclinações

novo_vetor_beta <- function(n,c,incl,b0=10,b1=1){
  x = 1:n
  x_c <- x[floor(c*n)]
  alpha1 = rCAT::rad2deg(atan(b1)) ; if(incl*alpha1 >= 90 |incl*alpha1 <= -90 ) stop("Inclinação muito acentuada -> novo angulo maior que 90º")
  alpha2 = rCAT::deg2rad(incl*alpha1) # Novo ângulo com efeito da inclinação desejada
  b1_new = tan(alpha2) # Valor do coeficiente angular da nova reta
  b2 = b1_new - b1
  b0_new = b0 - b2*x_c
  return(data.frame("Primeira_reta"=c(b0,b1),"Segunda_reta"=c(b0_new,b1_new)))
}

param <- expand.grid("n" = c(50 , 163 , 275 , 388 , 500),
                     "c" = c(0.25,0.5,0.75),
                     "incl"=c(0.05,0.10,0.15,0.25,0.50))

# ta errado a forma como os angulos são gerados....

tentar_corrigir <- function(n,c,b0=10,b1=1){
  
  
  vec_incl <- c(0.05,0.10,0.15,0.25,0.50) +1
  x = 1:n
  curve(b0+ b1*x,from=1,to=n)
  x_c <- x[floor(c*n)]
  incl <- vec_incl[1]
  alpha1 = rCAT::rad2deg(atan(b1)) ; if(incl*alpha1 >= 90 |incl*alpha1 <= -90 ) stop("Inclinação muito acentuada -> novo angulo maior que 90º")
  alpha2 = rCAT::deg2rad(incl*alpha1) # Novo ângulo com efeito da inclinação desejada
  b1_new = tan(alpha2) # Valor do coeficiente angular da nova reta
  b2 = b1_new - b1
  b0_new = b0 - b2*x_c
  retorno <- data.frame("Primeira_reta"=c(b0,b1),"Segunda_reta1"=c(b0_new,b1_new))
  curve(b0_new+ b1_new*x,from=1,to=n,col="yellow",add = TRUE)
  incl <- vec_incl[2]
  alpha1 = rCAT::rad2deg(atan(b1)) ; if(incl*alpha1 >= 90 |incl*alpha1 <= -90 ) stop("Inclinação muito acentuada -> novo angulo maior que 90º")
  alpha2 = rCAT::deg2rad(incl*alpha1) # Novo ângulo com efeito da inclinação desejada
  b1_new = tan(alpha2) # Valor do coeficiente angular da nova reta
  b2 = b1_new - b1
  b0_new = b0 - b2*x_c
  retorno <- cbind(retorno,"Segunda_reta2"=c(b0_new,b1_new))
  curve(b0_new+ b1_new*x,from=1,to=n,col="orange",add = TRUE)
  incl <- vec_incl[3]
  alpha1 = rCAT::rad2deg(atan(b1)) ; if(incl*alpha1 >= 90 |incl*alpha1 <= -90 ) stop("Inclinação muito acentuada -> novo angulo maior que 90º")
  alpha2 = rCAT::deg2rad(incl*alpha1) # Novo ângulo com efeito da inclinação desejada
  b1_new = tan(alpha2) # Valor do coeficiente angular da nova reta
  b2 = b1_new - b1
  b0_new = b0 - b2*x_c
  retorno <- cbind(retorno,"Segunda_reta3"=c(b0_new,b1_new))
  curve(b0_new+ b1_new*x,from=1,to=n,col="red",add = TRUE)
  incl <- vec_incl[4]
  alpha1 = rCAT::rad2deg(atan(b1)) ; if(incl*alpha1 >= 90 |incl*alpha1 <= -90 ) stop("Inclinação muito acentuada -> novo angulo maior que 90º")
  alpha2 = rCAT::deg2rad(incl*alpha1) # Novo ângulo com efeito da inclinação desejada
  b1_new = tan(alpha2) # Valor do coeficiente angular da nova reta
  b2 = b1_new - b1
  b0_new = b0 - b2*x_c
  retorno <- cbind(retorno,"Segunda_reta4"=c(b0_new,b1_new))
  curve(b0_new+ b1_new*x,from=1,to=n,col="purple",add = TRUE)
  incl <- vec_incl[5]
  alpha1 = rCAT::rad2deg(atan(b1)) ; if(incl*alpha1 >= 90 |incl*alpha1 <= -90 ) stop("Inclinação muito acentuada -> novo angulo maior que 90º")
  alpha2 = rCAT::deg2rad(incl*alpha1) # Novo ângulo com efeito da inclinação desejada
  b1_new = tan(alpha2) # Valor do coeficiente angular da nova reta
  b2 = b1_new - b1
  b0_new = b0 - b2*x_c
  retorno <- cbind(retorno,"Segunda_reta5"=c(b0_new,b1_new))
  curve(b0_new+ b1_new*x,from=1,to=n,col="blue",add = TRUE)
  
  return(retorno)
}
tentar_corrigir(50,0.5)
 

#  para o apendice --------------------------------------------------------


# c = 0.5
#tab_c50_A = cbind(Métrica = "A", new_results_epv1[which(new_results_epv1[,1] == "0.50"),-1]) 
#tab_c50_B = cbind(Métrica = "B", new_results_epv2[which(new_results_epv2[,1] == "0.50"),-1])
#print(xtable(rbind(tab_c50_A, tab_c50_B)), include.rownames=FALSE)
tab_c50_A = new_results_epv1[which(new_results_epv1[,1] == "0.5"),-1]
tab_c50_B = new_results_epv2[which(new_results_epv2[,1] == "0.5"),-1]

print(xtable(tab_c50_A), include.rownames=FALSE)
print(xtable(tab_c50_B), include.rownames=FALSE)

# c = 0.75
# tab_c75_A = cbind(Métrica = "A", new_results_epv1[which(new_results_epv1[,1] == "0.75"),-1]) 
# tab_c75_B = cbind(Métrica = "B", new_results_epv2[which(new_results_epv2[,1] == "0.75"),-1])
# print(xtable(rbind(tab_c75_A, tab_c75_B)), include.rownames=FALSE)

tab_c75_A = new_results_epv1[which(new_results_epv1[,1] == "0.75"),-1]
tab_c75_B = new_results_epv2[which(new_results_epv2[,1] == "0.75"),-1]

print(xtable(tab_c75_A), include.rownames=FALSE)
print(xtable(tab_c75_B), include.rownames=FALSE)

# c = 0.25

tab_c25_A = new_results_epv1[which(new_results_epv1[,1] == "0.25"),-1]
tab_c25_B = new_results_epv2[which(new_results_epv2[,1] == "0.25"),-1]

print(xtable(tab_c25_A), include.rownames=FALSE)
print(xtable(tab_c25_B), include.rownames=FALSE)




# Comparando métricas
# fixa n == 50
n50_metA = cbind(Métrica = "A", new_results_epv1[which(new_results_epv1[,2] == " 50"),-2]) 
n50_metB = cbind(Métrica = "B", new_results_epv2[which(new_results_epv2[,2] == " 50"),-2])
tbl_n50 = rbind(n50_metA, n50_metB)
r2 = tbl_n50[,3]
ct = tbl_n50[,2]
tbl_n50_v2 = tbl_n50[order(r2, ct),]
print(xtable(tbl_n50_v2), include.rownames=FALSE)

# fixa n == 163
n163_metA = cbind(Métrica = "A", new_results_epv1[which(new_results_epv1[,2] == "163"),-2])
n163_metB = cbind(Métrica = "B", new_results_epv2[which(new_results_epv2[,2] == "163"),-2])
tbl_n163 = rbind(n163_metA, n163_metB)
r2 = tbl_n163[,3]
ct = tbl_n163[,2]
tbl_n163_v2 = tbl_n163[order(r2, ct),]
print(xtable(tbl_n163_v2), include.rownames=FALSE)

# fixa n == 275
n275_metA = cbind(Métrica = "A", new_results_epv1[which(new_results_epv1[,2] == "275"),-2])
n275_metB = cbind(Métrica = "B", new_results_epv2[which(new_results_epv2[,2] == "275"),-2])
tbl_n275 = rbind(n275_metA, n275_metB)
r2 = tbl_n275[,3]
ct = tbl_n275[,2]
tbl_n275_v2 = tbl_n275[order(r2, ct),]
print(xtable(tbl_n275_v2), include.rownames=FALSE)

# fixa n == 388
n388_metA = cbind(Métrica = "A", new_results_epv1[which(new_results_epv1[,2] == "388"),-2])
n388_metB = cbind(Métrica = "B", new_results_epv2[which(new_results_epv2[,2] == "388"),-2])
tbl_n388 = rbind(n388_metA, n388_metB)
r2 = tbl_n388[,3]
ct = tbl_n388[,2]
tbl_n388_v2 = tbl_n388[order(r2, ct),]
print(xtable(tbl_n388_v2), include.rownames=FALSE)



# fixa n == 500
n500_metA = cbind(Métrica = "A", new_results_epv1[which(new_results_epv1[,2] == "500"),-2])
n500_metB = cbind(Métrica = "B", new_results_epv2[which(new_results_epv2[,2] == "500"),-2])
tbl_n500 = rbind(n500_metA, n500_metB)
r2 = tbl_n500[,3]
ct = tbl_n500[,2]
tbl_n500_v2 = tbl_n500[order(r2, ct),]
print(xtable(tbl_n500_v2), include.rownames=FALSE)

# para o texto ---------------------------------------------------------


# Tabelas de campeões dos cenários

# Métrica A
tbl_metA = resultado_epv[[1]]
#tbl_metA$value = as.numeric(as.character(tbl_metA$value))

summary_metA_50 = apply(tbl_metA[tbl_metA[,2]==50,], 1, function(x){
  valor_min = min(x[5:13])
  retorno = c(x[1:4], ifelse((x[5:13] == valor_min) == T, "X", ""))
  return(retorno)
})
summary_metA_163 = apply(tbl_metA[tbl_metA[,2]==163,], 1, function(x){
  valor_min = min(x[5:13])
  retorno = c(x[1:4], ifelse((x[5:13] == valor_min) == T, "X", ""))
  return(retorno)
})
summary_metA_275 = apply(tbl_metA[tbl_metA[,2]==275,], 1, function(x){
  valor_min = min(x[5:13])
  retorno = c(x[1:4], ifelse((x[5:13] == valor_min) == T, "X", ""))
  return(retorno)
})
summary_metA_388 = apply(tbl_metA[tbl_metA[,2]==388,], 1, function(x){
  valor_min = min(x[5:13])
  retorno = c(x[1:4], ifelse((x[5:13] == valor_min) == T, "X", ""))
  return(retorno)
})
summary_metA_500 = apply(tbl_metA[tbl_metA[,2]==500,], 1, function(x){
  valor_min = min(x[5:13])
  retorno = c(x[1:4], ifelse((x[5:13] == valor_min) == T, "X", ""))
  return(retorno)
})
#View(t(summary_metA))
print(xtable(t(summary_metA_50)[,-2]),include.rownames=FALSE)
print(xtable(t(summary_metA_163)[,-2]),include.rownames=FALSE)
print(xtable(t(summary_metA_275)[,-2]),include.rownames=FALSE)
print(xtable(t(summary_metA_388)[,-2]),include.rownames=FALSE)
print(xtable(t(summary_metA_500)[,-2]),include.rownames=FALSE)

# Métrica B
tbl_metB = resultado_epv[[2]]

summary_metB_50 = apply(tbl_metB[tbl_metB[,2]==50,], 1, function(x){
  valor_min = min(x[5:13])
  retorno = c(x[1:4], ifelse((x[5:13] == valor_min) == T, "X", ""))
  return(retorno)
})
summary_metB_163 = apply(tbl_metB[tbl_metB[,2]==163,], 1, function(x){
  valor_min = min(x[5:13])
  retorno = c(x[1:4], ifelse((x[5:13] == valor_min) == T, "X", ""))
  return(retorno)
})
summary_metB_275 = apply(tbl_metB[tbl_metB[,2]==275,], 1, function(x){
  valor_min = min(x[5:13])
  retorno = c(x[1:4], ifelse((x[5:13] == valor_min) == T, "X", ""))
  return(retorno)
})
summary_metB_388 = apply(tbl_metB[tbl_metB[,2]==388,], 1, function(x){
  valor_min = min(x[5:13])
  retorno = c(x[1:4], ifelse((x[5:13] == valor_min) == T, "X", ""))
  return(retorno)
})
summary_metB_500 = apply(tbl_metB[tbl_metB[,2]==500,], 1, function(x){
  valor_min = min(x[5:13])
  retorno = c(x[1:4], ifelse((x[5:13] == valor_min) == T, "X", ""))
  return(retorno)
})
#View(t(summary_metA))
print(xtable(t(summary_metB_50)[,-2]),include.rownames=FALSE)
print(xtable(t(summary_metB_163)[,-2]),include.rownames=FALSE)
print(xtable(t(summary_metB_275)[,-2]),include.rownames=FALSE)
print(xtable(t(summary_metB_388)[,-2]),include.rownames=FALSE)
print(xtable(t(summary_metB_500)[,-2]),include.rownames=FALSE)


# outro pro apendice ------------------------------------------------------

# Tabelas de campeões dos cenários

# Métrica A
tbl_metA = resultado_epv[[1]]
#tbl_metA$value = as.numeric(as.character(tbl_metA$value))

summary_metA_50 = tbl_metA[tbl_metA[,2]==50,]
summary_metA_163 = tbl_metA[tbl_metA[,2]==163,]
summary_metA_275 = tbl_metA[tbl_metA[,2]==275,]
summary_metA_388 = tbl_metA[tbl_metA[,2]==388,]
summary_metA_500 = tbl_metA[tbl_metA[,2]==500,]
#View(t(summary_metA))
print(xtable(t(summary_metA_50)[,-2]),include.rownames=FALSE)
print(xtable(t(summary_metA_163)[,-2]),include.rownames=FALSE)
print(xtable(t(summary_metA_275)[,-2]),include.rownames=FALSE)
print(xtable(t(summary_metA_388)[,-2]),include.rownames=FALSE)
print(xtable(t(summary_metA_500)[,-2]),include.rownames=FALSE)

# Métrica B
tbl_metB = resultado_epv[[2]]

summary_metB_50 = tbl_metB[tbl_metB[,2]==50,]
summary_metB_163 = tbl_metB[tbl_metB[,2]==163,]
summary_metB_275 = tbl_metB[tbl_metB[,2]==275,]
summary_metB_388 = tbl_metB[tbl_metB[,2]==388,]
summary_metB_500 = tbl_metB[tbl_metB[,2]==500,]
#View(t(summary_metA))
print(xtable(t(summary_metB_50)[,-2]),include.rownames=FALSE)
print(xtable(t(summary_metB_163)[,-2]),include.rownames=FALSE)
print(xtable(t(summary_metB_275)[,-2]),include.rownames=FALSE)
print(xtable(t(summary_metB_388)[,-2]),include.rownames=FALSE)
print(xtable(t(summary_metB_500)[,-2]),include.rownames=FALSE)



# 8.2) Nova Visualização Resultados c < 1 ------------------------------------------

#setwd("C:/Users/Igor/OneDrive/Documentos/TCC/Códigos Adicionais")
library(reshape2)
library(lattice)
library(data.table)

n = c(50 , 163 , 275 , 388 , 500)
c = c(0.25,0.5,0.75)
r_quad = c(0.2,0.5,0.7,0.85,0.95)
met = c("A","B")

# métrica A
resultado_epvA <- metodoA
resultado_epvB <- metodoB
resultado_epv =list(resultado_epvA, resultado_epvB)


for (k in 1:2){
  
  for (ci in 1:3){
    
    for (r2 in 1:5){
      
      for (ni in 1:5){
        
        resultado_filtro <-  as.data.frame(resultado_epv[[k]][which(resultado_epv[[k]][,1] == c[ci] & resultado_epv[[k]][,2] == n[ni] & resultado_epv[[k]][,3] == r_quad[r2]),])
        
        resultado_filtro$Rec.CUSUM <- as.vector(as.list(resultado_filtro$Rec.CUSUM))
        resultado_filtro$OLS.CUSUM <- as.vector(as.list(resultado_filtro$OLS.CUSUM))
        resultado_filtro$fluctuation <- as.vector(as.list(resultado_filtro$fluctuation))
        resultado_filtro$Rec.MOSUM <- as.vector(as.list(resultado_filtro$Rec.MOSUM))
        resultado_filtro$OLS.MOSUM <- as.vector(as.list(resultado_filtro$OLS.MOSUM))
        resultado_filtro$ME <- as.vector(as.list(resultado_filtro$ME))
        resultado_filtro$supF <- as.vector(as.list(resultado_filtro$supF))
        resultado_filtro$aveF <- as.vector(as.list(resultado_filtro$aveF))
        resultado_filtro$expF <- as.vector(as.list(resultado_filtro$expF))
        setDT(resultado_filtro)
        
        result_melt <-  melt(resultado_filtro[,-(1:3)], id = "incl")
        result_melt$value = as.numeric(result_melt$value)
        result_melt$incl = as.factor(as.numeric(result_melt$incl))
        
        plot.save <-  xyplot(value ~ incl,
                             data = result_melt,
                             groups = variable,
                             type = 'l',
                             lwd = 2,
                             auto.key = list(title = "Testes", lwd=4, cex=1, space = 'right'),
                             main = paste("Métrica:", met[k], "\n", "c =", c[ci], "\n", "n =", n[ni], "\n", "r2 =", r_quad[r2])
        )
        trellis.device(device="png", filename=paste("xyplot_corte",c[ci]*100,"_","metrica",met[k],"_","n",n[ni],"_","rquad",100*r_quad[r2],".png", sep = ""))
        print(plot.save)
        dev.off()
      }
    }
  }
}











