### Projeto: nimcno
### Especificacao das condicoes iniciais

#tab.stationary(diff(tsm[,c("logpibr", "logm1")]))


## Epecifica modelos individuais

# estima cada modelo separadamente

library(nimcno)

bn.ipca <- ucmodel(x = macro95[,"ipca"], l = 3,
                   init = c(rep(c(1.5, 0.4, 0.2), 1),# ar pars
                            c(-0.5), # var arima
                            c(-1), # var level
                            rep(0.01, 1) # cov arima level
                            ))

ucplot(bn.ipca, state = "level")
core <- bn.ipca$out$alphahat[,'level']

bn.ipca$fit$model$T
bn.ipca$fit$model$Q


# selicr
x <- tsm[,"selicr"]
bn.selicr <- ucmodel(x = x, l = 3,
                   init = c(rep(c(1.5, 0.4, 0.2), 1),# ar pars
                            c(-0.3), # var arima
                            c(-1.6), # var level
                            rep(0.01, 1) # cov arima level
                   ))

tsplot(bn.selicr$out$alphahat[,'level'], tsm[,"selicr"])
bn.selicr$fit$model$T
bn.selicr$fit$model$Q

# logpibr
bn.logpibr <- ucmodel(x = (macro95[,"logpibr"]), l = 2,
                      init = c(rep(c(2, 0.2), 1),# ar pars
                               c(-2), # var arima
                               c(-3), # var level
                               rep(0.01, 1) # cov arima level
                      ))

ucplot(bn.logpibr, state = "level")
tsplot(bn.logpibr$out$alphahat[,'arima1'])
#head(bn.logpibr$out$alphahat)
bn.logpibr$fit$model$T
bn.logpibr$fit$model$Q


# logm1
bn.logm1 <- ucmodel(x = tsm[,c('logm1')], l=3, 
                    init=c(1.5, -0.6, 0.02, 
                           -0.5,  # var arima
                           -1.7, # var level
                           0.01 # cov arima level
                           ))
tsplot(bn.logm1$out$alphahat[,'level'], tsm[,"logm1"])
tsplot(bn.logm1$out$alphahat[,'arima1'])
bn.logm1$fit$model$T
bn.logm1$fit$model$Q

exp(2*-1.7)/exp(2*-0.5)

## Especifica o modelo conjunto sem correlacao entre arima e level

x <- macro95[,c('ipca', 'selicr', 'logpibr', 'logm1')]
bn1 <- ucmodel(x = x, l=3, 
              init=c(c(0.70, -0.19, 0.06),
                     c(0.71, -0.11, 0.08),
                     c(0.63, -0.23, 0.08),
                     c(0.73, 0.72, -0.99),
                     c(-0.5, -0.3, -2, -1),  # var arima
                     c(-1, -1.6, -3, -2), # var level
                     rep(0.01, 6), # cov arima 
                     rep(0.01, 6) # cov level
                     #rep(0.01, 16) # cov arima level
                     ), corre = FALSE)

y <- bn1$out$alphahat
y <- y[,grepl('level', colnames(y))]
tsplot(y, x)

x <- macro95[,c('ipca', 'selicr', 'logpibr', 'logm1')]
bn2 <- ucmodel(x = x, l=3, 
               init=c(c(0.47, 0.05, -0.04),
                      c(0.61, 0.01, -0.03),
                      c(0.29, -0.25, 0.14),
                      c(1.06, -0.07, -0.40),
                      c(-0.5, -0.3, -1.5, -1),  # var arima
                      c(-1, -1.3, -3.3, -2.3), # var level
                      rep(0.01, 6), # cov arima 
                      rep(0.01, 6), # cov level
                      rep(0.01, 16) # cov arima level
               ), corre = TRUE)

y <- bn2$out$alphahat
y <- y[,grepl('level', colnames(y))]
tsplot(y, x)



x = macro95[,c('ipca', 'selicr')]
library(nimcno)
bn2 <- ucmodel(x = macro95[,c('ipca', 'selicr')], l=2, 
               init=c(c(2, 0.2),
                      c(2, 0.2),
                      c(-0.5, -0.5),  # var arima
                      c(-1.2, -1.2), # var level
                      0.001, # cov level 
                      0.01,  # cov arima
                      rep(0.001, 4) # cov level arima
               ), corre = TRUE)

ucplot(bn2, state = "level", ic = 0.95)
bn2$fit$model$T
ldl(bn2$fit$model$Q[,,1])

-0.8684986

### Calcula o intervalo de confiança 

teste <- ic.ucmodel(bn2, state = "level", ic = 0.95) 
teste$upper


x <- macro95[,c('ipca', 'selicr', 'logpibr')]
bn2 <- ucmodel(x = x, l=3, 
               init=c(c(0.47, 0.05, -0.04),
                      c(0.61, 0.01, -0.03),
                      c(0.29, -0.25, 0.14),
                      c(1.06, -0.07, -0.40),
                      c(-0.5, -0.3, -1.5, -1),  # var arima
                      c(-1, -1.3, -3.3, -2.3), # var level
                      rep(0.01, 6), # cov arima 
                      rep(0.01, 6), # cov level
                      rep(0.01, 16) # cov arima level
               ), corre = TRUE)


bn2 <- ucmodel(x = macro95[,c('ipca', 'selicr', 'logpibr')], l=2, 
               init=c(c(2, 0.2),
                      c(2, 0.2),
                      c(2, 0.2),
                      c(-0.5, -0.5, -2),  # var arima
                      c(-1.3, -1.3, -3), # var level
                      c(0.1, 0.1, 0.1), # cov level 
                      c(0.1, 0.01, 0.01),  # cov arima
                       c(0.001, 0.001, 0.001),# cov level arima
                       c(0.001, 0.001, 0.001),
                       c(0.001, 0.001, 0.001)
               ), corre = TRUE)

ucplot(bn2, state = "level", ic = 0.95)
bn2$fit$model$T
KFAS::ldl(bn2$fit$model$Q[,,1])
pars.ucmodel(bn2)

ucplot(bn2, state = "slope")
tsplot(bn2$out$alphahat[,"slope.selicr"])


x <- cbind(bn.ipca$out$alphahat[,"level"],
           bn2$out$alphahat[,"level.ipca"],
           bn2$out$alphahat[,"level.selicr"])
y <- cbind(zoo::as.zoo(bn.ipca$out$model$y),
           zoo::as.zoo(bn2$out$model$y))
dp <- cbind(ic.ucmodel(bn.ipca, state = "level")$v,
            ic.ucmodel(bn2, state = "level")$v)
ic <- ic(x, dp = dp, nc = nc)

colnames(x) <- c("IPCA univariado", "IPCA multivariado", "Selicr multivariado")

