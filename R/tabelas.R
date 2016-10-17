#########
### Projeto: nimcno
########

#' Significancia de um teste
#' 
#' Fornece o valor da estatistica de um teste junto com a indicacao
#' da significancia para valores 10\%, 5\% e 1\%.
#' 
#' @param stat valor da estatistica do teste
#' @param cval vetor com os valores criticos do teste para 10\%, 5\% e 1\%
#' 
#' @return um objeto "character" com valor da estatistica do teste formatado
#' junto com "***" para 1\%, "**" para 5\%, "*" para 10\% e " " nao significativo
#' 
#' @export
#'

sig <- function(stat, cval){
  x <- sum(abs(stat)>abs(cval))
  if(x==3) sig <- "***"
  if(x==2) sig <- "** "
  if(x==1) sig <- "*  "
  if(x==0) sig <- "   "
  return(c(paste(format(round(stat,3), digits = 3, nsmall = 3, decimal.mark=","), sig, sep = "")))
}

#' Tabela para o teste de estacionaridade ADF e KPSS
#' 
#' Efetua o teste ADF e KPSS para series em nivel e 
#' diferenciadas e fornece uma tabela com as estatisticas
#' 
#' @param x um objeto mts com series temporais multivariadas
#' @param d logico, TRUE (padrao) indica que teste deve ser 
#' efetuado tambem para series estacionarias
#' 
#' @return tabela tabela com as estatisticas dos testes
#' 
#' @export
#' 


tab.stationary <- function(x, d = TRUE){
  n <- dim(x)[2]
  tab <- data.frame(Variavel=colnames(x), tendencia="sim", ADF.lag=NA, ADF=NA, KPSS=NA)
  for(i in 1:n){
    tmp <- try(stats::na.fail(x[,i]), silent = TRUE)
    if(inherits(tmp, "try-error")){
      y <- stats::na.omit(x[,i])
      warning("NAs omited")
    }else{
      y <- x[,i]
    }
    k <- adf.lag(y, type = "trend", lag.max = 15, selectlags = "AIC")$lag
    adf <- urca::ur.df(y, type = "trend", lags=k)
    kpss <- urca::ur.kpss(y, type = "tau")
    tab[i,'ADF.lag'] <- k
    tab[i,'ADF'] <- sig(stat = adf@teststat[1], cval = adf@cval[1,])
    tab[i,'KPSS'] <- sig(stat = kpss@teststat, cval = kpss@cval[-3])
  }
  if(d){
    tab2 <- data.frame(Variavel=paste0("diff.",colnames(x)), tendencia="nao", ADF.lag=NA, ADF=NA, KPSS=NA)
    for(i in 1:n){
      y <- na.omit(diff(x[,i]))
      k <- adf.lag(y, type = "drift", lag.max = 15, selectlags = "AIC")$lag
      adf <- urca::ur.df(y, type = "drift", lags=k)
      kpss <- urca::ur.kpss(y, type = "mu")
      tab2[i,'ADF.lag'] <- k
      tab2[i,'ADF'] <- sig(stat = adf@teststat[1], cval = adf@cval[1,])
      tab2[i,'KPSS'] <- sig(stat = kpss@teststat, cval = kpss@cval[-3])
    }
    tabela <- rbind(tab,tab2)
    return(tabela)
  }
  tabela <- tab
  return(tabela)
}

#' Teste de cointegracao de Johansen
#' 
#' Conducts the Johansen procedure on a given data set. 
#' The "trace" or "eigen" statistics are reported. 
#' A wrapper for the \link[urca]{ca.jo} function in the urca package. 
#' 
#' 
#' @param x Data matrix to be investigated for cointegration
#' @param ecdet Character, "none" for no intercept in 
#' cointegration, "const" for constant term in cointegration 
#' and "trend" for trend variable in cointegration.
#' @param K	The lag order of the series (levels) in the VAR
#' @param season	If seasonal dummies should be included, the
#'  data frequency must be set accordingly, i.e "4" for 
#'  quarterly data.
#' 
#' @return A list with the following elements:
#'   \item{\code{eigen}}{estatistica do teste do lambda-max}
#'   \item{\code{trace}}{estatistica do teste do traco}
#'   
#'   @import urca
#'   
#'   @export
#'   


cointeg <- function(x, ecdet, K, season = NULL){
  eigen <- urca::ca.jo(x, type = "eigen", ecdet = ecdet, K = K, season = season)
  trace <- urca::ca.jo(x, type = "trace", ecdet = ecdet, K = K, season = season)
  
  tmp <- cbind(eigen@teststat,eigen@cval)
  sig.eigen <- apply(tmp, 1, function(x) sig(stat = x[1], cval = x[-1]))
  tmp <- cbind(trace@teststat, trace@cval)
  sig.trace <- apply(tmp, 1, function(x) sig(stat = x[1], cval = x[-1]))
  return(list(eigen=sig.eigen, trace=sig.trace))
}


#' Tabela para o teste de cointegracao 
#' 
#' Efetua o teste de cointegracao de Johansen para series e 
#' fornece uma tabela com as estatisticas 
#' 
#' @param x um objeto mts com series temporais multivariadas
#' 
#' @return tabela com as estatisticas dos testes
#' 
#' @export
#' 

tab.cointeg <- function(x){
  y <- cointeg(x, ecdet = "const", K = 4)
  n <- dim(x)[2]
  tab <- data.frame(Posto = 0:(n-1), Traco = y$trace[n:1],
                    Lambda.max = y$eigen[n:1])
  return(tab)
}





#' Testa os criterios de Marques et al (2003)
#'
#' Testa os criterios de Marques et al (2003) que uma 
#' medida de nucleo da inflacao deve atender
#'
#' @param x serie do nucleo
#' @param y serie da inflacao 
#' 
#' @return  um vetor com os teste "ADF", "t alpha", "t gamma", 't lambda', 'F thetas'
#'
#' @import zoo dynlm urca lmtest
#' 
#' @export
#' 


marques <- function(x,y){
  requireNamespace("zoo")
  requireNamespace("dynlm")
  
  if(sum(round(stats::tsp(x),  3)!=round(stats::tsp(y), 3))==1){
    warning("series com inicio, fim ou frequencia diferentes, 
            usando somente a cobertura temporal comum")
    ini <- max(stats::tsp(x)[1], stats::tsp(y)[1])
    fim <- min(stats::tsp(x)[2], stats::tsp(y)[2])
    x <- stats::window(x, start = ini, end = fim)
    y <- stats::window(y, start = ini, end = fim)
  }
  y <- zoo::as.zoo(y)
  x <- zoo::as.zoo(x)
  
  ## Teste adf 
  z <- y - x
  fit.adf <- adf.lag(z, type = "drift", lag.max = 15, selectlags = "AIC")
  k <- fit.adf$lag
  adf <- urca::ur.df(z, lags=k, type= "drift")
  t.adf <- sig(stat = adf@teststat[1], cval = adf@cval[1,])
  t.alpha <- lmtest::coeftest(fit.adf$fit)[1,4]
  
  ## estima os mecanismos de correcao de erro para a inflacao y e
  ## e o nucleo x
  
  # escolha as defasagens k por AIC
  result.y <- result.x <- vector()
  for(i in 1:15){
    fity <- dynlm::dynlm(d(y) ~ L(I(y-x), 1) + L(d(y), 1:i) + L(d(x), 1:i))
    fitx <- dynlm::dynlm(d(x) ~ L(I(y-x), 1) + L(d(y), 1:i) + L(d(x), 1:i))
    result.y[i] <- stats::AIC(fity)
    result.x[i] <- stats::AIC(fitx)
  }
  k <- which.min(result.y)
  fity <- dynlm::dynlm(d(y) ~ L(I(y-x), 1) + L(d(y), 1:k) + L(d(x), 1:k))
  k <- which.min(result.x)
  fitx <- dynlm::dynlm(d(x) ~ L(I(y-x), 1) + L(d(y), 1:k) + L(d(x), 1:k))
  
  # exogeneidade fraca
  # p valor do teste t sobre gamma (y) e lambda (x)
  #t.gamma.y <- lmtest::coeftest(fity, vcov. = sandwich::NeweyWest(fity))[2,4]
  t.gamma.y <- summary(fity)$coefficients[2,4]
  #t.lambda.x <- lmtest::coeftest(fitx, vcov. = sandwich::NeweyWest(fitx))[2,4]
  t.lambda.x <- summary(fitx)$coefficients[2,4]
  
  # exogeneidade forte
  # p valor teste F sobre os thethas
  k <- which.min(result.x)
  fit1 <- dynlm::dynlm(d(x) ~ L(d(y), 1:k) + L(d(x), 1:k))
  fit2 <- dynlm::dynlm(d(x) ~ L(d(x), 1:k))
  #test.f <- lmtest::waldtest(fit1, fit2, test = "F", vcov = sandwich::NeweyWest(fit1))
  test.f <- lmtest::waldtest(fit1, fit2, test = "F")
  f.theta <- test.f$`Pr(>F)`[2]
  
  pvalue <- c(t.alpha, t.gamma.y, t.lambda.x, f.theta)
  result <- c(t.adf, n2tab(pvalue))
  names(result) <- c("ADF", "t alpha", "t gamma y", 't lambda x', 'F thetas')
  return(result)
}

#' Tabela para as condicoes Marques et. al. (2003)
#' 
#' Efetua o teste para as condicoes Marques et. al. (2003) para o nucleo e inflacao
#'  e fornece uma tabela com as estatisticas 
#' 
#' @param y um objeto ts com a inflacao observada
#' @param x um objeto mts com series dos nucleos a serem testados
#' 
#' @return tabela com as estatisticas dos testes
#' 
#' @export
#' 

tab.marques <- function(y, x){
  n <- dim(x)[2]
  tab <- data.frame(nucleos=colnames(x), ADF=NA, t.alpha=NA, t.gamma=NA,
                    t.lambda=NA, F.thetas=NA)
  for(i in 1:n){
    tab[i,2:6] <- marques(y=y,x=stats::na.omit(x[,i]))
  }
  return(tab)
}

#' Transforma numero e texto formatado
#' 
#' Transforma numero e texto formatado para serem mostrados na
#' tabela, o numero arredondado para 3 casas decimais e separado
#' por virgula
#' 
#' @param x um objeto com valores numericos
#' 
#' @return o mesmo objeto que foi fornecido com numeros 
#' transformados em texto formatado
#' 
#' @export
#' 

n2tab <- function(x){
  format(round(x, 3), digits = 3, nsmall = 3, decimal.mark=",") 
}




#' Seleciona a defasagem para o teste ADF
#' 
#' Seleciona a defasagem para o teste ADF segundo o criterio
#' AIC ou BIC para um maximo de defasagens
#' 
#' @param y The vector tested for a unit root
#' @param type Test type, either "none", "drift" or "trend".
#' @param lag.max The maximum number of lags considered
#' @param selectlags Lag selection can be achieved according to the Akaike "AIC" or the Bayes "BIC" information criteria. 
#' 
#' @return A defasagem selecionda
#' 
#' @import stats
#' 
#' @export

adf.lag <- function (y, type = c("none", "drift", "trend"), lag.max = 15, selectlags = c("AIC", "BIC")) 
{
  selectlags <- match.arg(selectlags)
  type <- match.arg(type)
  if (ncol(as.matrix(y)) > 1) 
    stop("\ny is not a vector or univariate time series.\n")
  if (any(is.na(y))) 
    stop("\nNAs in y.\n")
  y <- as.vector(y)
  lag <- as.integer(lag.max)
  if (lag < 0) 
    stop("\nLags must be set to an non negative integer value.\n")
  CALL <- match.call()
  DNAME <- deparse(substitute(y))
  x.name <- deparse(substitute(y))
  z <- diff(y)
  n <- length(z)
  x <- stats::embed(z, lag.max)
  z.diff <- x[, 1]
  z.lag.1 <- y[lag.max:n]
  tt <- lag.max:n
  if (lag.max > 1) {
    critRes <- rep(NA, lag.max)
    for (i in 2:(lag.max)) {
      z.diff.lag = x[, 2:i]
      if (type == "none") 
        result <- stats::lm(z.diff ~ z.lag.1 - 1 + z.diff.lag)
      if (type == "drift") 
        result <- stats::lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)
      if (type == "trend") 
        result <- stats::lm(z.diff ~ z.lag.1 + 1 + tt + z.diff.lag)
      critRes[i] <- stats::AIC(result, k = switch(selectlags, 
                                           AIC = 2, BIC = log(length(z.diff))))
    }
    lag <- which.min(critRes)
  }else{
    lag <- lag.max
  }
  z.diff.lag = x[, 2:lag]
  if (type == "none") 
    result <- stats::lm(z.diff ~ z.lag.1 - 1 + z.diff.lag)
  if (type == "drift") 
    result <- stats::lm(z.diff ~ z.lag.1 + 1 + z.diff.lag)
  if (type == "trend") 
    result <- stats::lm(z.diff ~ z.lag.1 + 1 + tt + z.diff.lag)
  return(list(lag = lag, fit = result))
}