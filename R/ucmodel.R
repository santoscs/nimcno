## Estima tendencia e ciclo com modelo 
## de componentes nao observados
##


#' Estima o modelo de componentes nao observaveis com tendencia
#' e ciclo por maxima verossimelhanca
#' 
#' Estima o modelo de componentes nao observaveis com tendencia
#' e ciclo por maxima verossimelhanca usando o pacote "KFAS"
#' 
#' @param x um objeto mts com a serie temporia multivariada
#' @param lag numero de defagens do modelo AR 
#' @param init valores iniciais dos parametros a serem estimados 
#' @param corre logical, se TRUE (padrao) modelo com componentes 
#' correlacionado, se FALSE modelo com componentes restrito se 
#' correlacao
#'
#' @return uma lista com os seguinte objetos fit (\link[KFAS]{fitSSM}) e out (\link[KFAS]{KFS})
#' 
#' @import KFAS
#' 
#' @export
#'       

ucmodel <- function(x, lag, init, corre = TRUE){
  requireNamespace("KFAS")
  p <- dim(as.matrix(x))[2]
  ar0 <- lapply(apply(matrix(0, lag, p), 2, list), unlist)
  model <- SSModel(x ~ -1 + SSMtrend(degree = 2, Q=list(diag(1,p), diag(0,p)), type = "distinct", index = 1:p) +
                     SSMarima(ar = ar0, Q=diag(1,p), index =1:p),
                   H = diag(0, p))
  likfn <- function(pars, model){
    arl <- lapply(apply(matrix(pars[1:(lag*p)], lag, p), 2, list),unlist)
    tmp <- try(SSModel(x ~ -1 + SSMtrend(degree = 2, Q=list(diag(exp(2*pars[((lag+1)*p+1):((lag+2)*p)]), p), diag(0,p)), type = "distinct") +
                         SSMarima(ar = lapply(arl, artransform), Q=diag(exp(2*pars[(lag*p+1):((lag+1)*p)]), p), index = 1:p),
                       H = diag(0, p)), silent = TRUE)
    if(!inherits(tmp, "try-error")){
      model["T"] <- tmp$T
      model["R"] <- tmp$R
      model["P1"] <- tmp$P1
      
      Q <- tmp$Q[,,1]
      types <- attributes(tmp)$eta_types
      #covariancia level
      snn <- pars[((lag+2)*p + 1):((lag+2)*p + ((p^2-p)/2))] #* sn[upper.tri(sn)]
      id1 <- types == "level"
      Q[id1,id1][upper.tri(Q[id1,id1])] <- snn
      # covariancia arima
      see <- pars[((lag+2)*p + 1 + ((p^2-p)/2)):((lag+2)*p + 2*((p^2-p)/2))] #*se[upper.tri(se)]
      id2 <- types == "arima"
      Q[id2,id2][upper.tri(Q[id2,id2])] <- see
      # covariancia level arima 
      if(corre){
        snene <- pars[((lag+2)*p + 1 + 2*((p^2-p)/2)):((lag+2)*p + 2*((p^2-p)/2) + p*p)] 
        Q[id1, id2] <- snene
      }
      model["Q"][ , , 1] <- crossprod(Q)
      model
    } else {
      model
    }
  }
  fit <- fitSSM(model, inits = init,
                updatefn = likfn, 
                method = "BFGS",
                hessian = TRUE,
                maxiter = 5000,
                control = list(reltol = 1e-07))
  
  if(fit$optim.out$convergence!=0) warning("A otimizacao da funcao verossimilhanca nao convergiu")
  out <- KFS(fit$model)
  return(list(out=out, fit=fit))
}




#' Cria o modelo de componentes nao observaveis com tendencia
#' e ciclo 
#' 
#' Cria o modelo de componentes nao observaveis com tendencia
#' usando o pacote "KFAS"
#' 
#' @param x um objeto mts com a serie temporia multivariada
#' @param lag numero de defagens do modelo AR 
#' @param pars valores dos parametros  
#' @param corre logical, se TRUE (padrao) modelo com componentes 
#' correlacionado, se FALSE modelo com componentes restrito se 
#' correlacao
#'
#' @return uma lista com os seguinte objetos fit (\link[KFAS]{fitSSM}) e out (\link[KFAS]{KFS})
#' 
#' @import KFAS
#' 
#' @export
#'       

ucm <- function(x, lag, pars, corre = TRUE){
  requireNamespace("KFAS")
  p <- dim(as.matrix(x))[2]
  arl <- lapply(apply(matrix(pars[1:(lag*p)], lag, p), 2, list),unlist)
  tmp <- try(SSModel(x ~ -1 + SSMtrend(degree = 2, Q=list(diag(exp(2*pars[((lag+1)*p+1):((lag+2)*p)]), p), diag(0,p)), type = "distinct") +
                       SSMarima(ar = lapply(arl, artransform), Q=diag(exp(2*pars[(lag*p+1):((lag+1)*p)]), p), index = 1:p),
                     H = diag(0, p)), silent = TRUE)
  if(!inherits(tmp, "try-error")){
    Q <- tmp$Q[,,1]
    types <- attributes(tmp)$eta_types
    #covariancia level
    snn <- pars[((lag+2)*p + 1):((lag+2)*p + ((p^2-p)/2))] #* sn[upper.tri(sn)]
    id1 <- types == "level"
    Q[id1,id1][upper.tri(Q[id1,id1])] <- snn
    # covariancia arima
    see <- pars[((lag+2)*p + 1 + ((p^2-p)/2)):((lag+2)*p + 2*((p^2-p)/2))] #*se[upper.tri(se)]
    id2 <- types == "arima"
    Q[id2,id2][upper.tri(Q[id2,id2])] <- see
    # covariancia level arima 
    if(corre){
      snene <- pars[((lag+2)*p + 1 + 2*((p^2-p)/2)):((lag+2)*p + 2*((p^2-p)/2) + p*p)] 
      Q[id1, id2] <- snene
    }
    tmp["Q"][ , , 1] <- crossprod(Q)
  } else {
    stop("parametros nao factiveis")
  }
  return(tmp)
}





#' Create a State Space Model Object Arima of Class SSModel
#' 
#' Adaptadacao da funcao \link[KFAS]{SSMarima} do package KFAS 
#' para permitir parametros regressivos diferentes entre as series
#' 
#' @param ar a list of numeric vector containing the autoregressive coeffients.
#' @param ma nao usado (sempre NULL).
#' @param d a degree of differencing (nao usado, sempre 0)
#' @param Q a p x p covariance matrix of the disturbances 
#' (or in the time varying case p x p x n array), where 
#' $p$ = length(index)
#' @param stationary logical value indicating whether a stationarity of the 
#' arima part is assumed. Defaults to TRUE.
#' @param index A vector indicating for which series the corresponding 
#' components are constructed.
#' @param n Length of the series, only used internally for dimensionality check.
#' @param ynames names of the times series, used internally.
#' 
#' @return Object of class SSModel (ver \link[KFAS]{SSModel})
#' 
#' @export
#' 

SSMarima <-  function (ar = NULL, ma = NULL, d = 0, Q, stationary = TRUE, 
                       index, n = 1, ynames) 
{ 
  if (!sapply(ar, is.null) && stationary && !sapply(ar, function (x) all(Mod(polyroot(c(1,-x))) > 1))) 
    stop("ARIMA part is non-stationary.")
  if (missing(index)) 
    index <- 1
  p <- length(index)
  if (!missing(ynames) && !is.null(ynames)) {
    ynames <- paste0(".", ynames)
  }
  else ynames <- ""
  if (missing(Q)) {
    Q <- diag(p)
  }
  else {
    if (length(Q) == 1) 
      Q <- matrix(Q)
    if (any(dim(Q)[1:2] != p) || length(dim(Q)) > 2) 
      stop("Misspecified Q, argument Q must be (p x p) matrix where p is the number of series.")
  }
  ar_length <- length(ar[[1]])
  ma_length <- length(ma)
  d <- max(d, 0)
  m1 <- max(ar_length, ma_length + 1) + d
  k <- p
  Z_univariate <- matrix(0, 1, m1)
  P1inf_univariate <- matrix(0, m1, m1)
  T_univariate2 <- vector('list',length(ar))
  for(i in 1:length(ar)){
    T_univariate2[[i]] <- matrix(0, m1, m1)
  }
  R_univariate <- matrix(0, m1, 1)
  Z_univariate[1, 1:(d + 1)] <- 1
  if (d > 0) {
    for(i in 1:length(ar)){
      T_univariate2[[i]][1:d, 1:d][upper.tri(T_univariate2[[i]][1:d, 1:d], 
                                             diag = TRUE)] <- 1
      T_univariate2[[i]][1:d, (d + 1)] <- 1
    }
    P1inf_univariate[1:d, 1:d] <- diag(1, d)
  }
  if (ar_length > 0)
    for(i in 1:length(ar)){
      T_univariate2[[i]][(d + 1):(d + ar_length), d + 1] <- ar[[i]]
    }
  if (m1 > (d + 1)) 
    for(i in 1:length(ar)){
      T_univariate2[[i]][(d + 1):(m1 - 1), (d + 2):m1] <- diag(1, max(ar_length, ma_length + 1) - 1)
    }
  R_univariate[d + 1, 1] <- 1
  if (ma_length > 0) 
    R_univariate[(d + 2):(d + 1 + ma_length)] <- ma
  m <- p * m1
  Z <- matrix(0, p, m)
  T <- P1 <- P1inf <- matrix(0, m, m)
  R <- matrix(0, m, p)
  for (i in 1:p) {
    Z[i, ((i - 1) * m1 + 1):(i * m1)] <- Z_univariate
    T[((i - 1) * m1 + 1):(i * m1), ((i - 1) * m1 + 1):(i * 
                                                         m1)] <- T_univariate2[[i]]
    R[((i - 1) * m1 + 1):(i * m1), i] <- R_univariate
    P1inf[((i - 1) * m1 + 1):(i * m1), ((i - 1) * m1 + 1):(i * 
                                                             m1)] <- P1inf_univariate
  }
  if (stationary) {
    nd <- which(diag(P1inf) == 0)
    mnd <- length(nd)
    temp <- try(solve(a = diag(mnd^2) - matrix(kronecker(T[nd, 
                                                           nd], T[nd, nd]), mnd^2, mnd^2), b = c(R[nd, , drop = FALSE] %*% 
                                                                                                   Q %*% t(R[nd, , drop = FALSE]))), TRUE)
    if (class(temp) == "try-error") {
      stop("ARIMA part is numerically too close to non-stationarity.")
    }
    else P1[nd, nd] <- temp
  }
  else diag(P1inf) <- 1
  state_names <- paste0(rep(paste0("arima", 1:m1), p), rep(ynames, 
                                                           each = m1))
  list(index = index, m = m, k = k, Z = Z, T = T, R = R, Q = Q, 
       a1 = matrix(0, m, 1), P1 = P1, P1inf = P1inf, tvq = 0, 
       tvr = 0, tvz = 0, state_names = state_names)
}



#' Extrair os paramentros de objetos ucmodel
#' 
#' Extrair os paramentros do modelo de componentes nao observaveis
#' estimado por ucmodel
#' 
#' @param x um objeto "ucmodel"
#' 
#' @return uma lista com ar (parametros autoregressivos), sigma
#' (vaiancias) e rho (covariancias)
#' 
#' @import utils
#' 
#' @export
#' 

pars.ucmodel <- function(x){
  # parametros autoregressivos
  ssT <- x$out$model$T[,,1]
  nome <- colnames(x$fit$model$y)
  types <- colnames(ssT)
  arpars <- NULL
  if(!is.null(nome)){
    for(i in 1:length(nome)){
      id1 <- as.logical(grepl("arima", types) * grepl(nome[i], types))
      id2 <- grepl(paste0("arima1.", nome[i]), types)
      arpars <- cbind(arpars, ssT[id1, id2])
    }
  } else{
    id1 <- grepl("arima", types)
    id2 <- grepl("arima1", types)
    arpars <- ssT[id1, id2]
  }
  
  # tipos de variancia
  ssQ <- x$out$model$Q[,,1]
  types <- attributes(x$out$model)$eta_types
  
  if(!is.null(nome)){
    # acresenta nome das series
    types[grep("level", types)] <- paste0(types[grep("level", types)], ".", nome)
    types[grep("slope", types)] <- paste0(types[grep("slope", types)], ".", nome)
    types[grep("arima", types)] <- paste0(types[grep("arima", types)], ".", nome)
  }
  colnames(ssQ) <- types
  rownames(ssQ) <- types
  
  
  # separar variancia e covariancia
  # level
  id1 <- grepl("level", types)
  cov.level <- ssQ[id1, id1]
  
  # arima
  id2 <- grepl("arima", types)
  cov.arima <- ssQ[id2, id2]
  
  # cruzada level arima
  cov.level.arima <- ssQ[id1, id2]
  
  
  #correlacao
  id <- grepl("arima", types) | grepl("level", types)
  cor.m <- cov2cor(ssQ[id,id])
  
  id1 <- grepl("level", colnames(cor.m))
  tab1 <- n2tab(cor.m[id1, id1])
  tab1[lower.tri(tab1)] <- ifelse(tab1[lower.tri(tab1)]>0, "+","-")
  id2 <- grepl("arima", colnames(cor.m))
  tab2 <- n2tab(cor.m[id2, id2])
  tab2[lower.tri(tab2)] <- ifelse(tab2[lower.tri(tab2)]>0, "+","-")
  tab3 <- n2tab(cor.m[id1, id2])
  tab4 <- ifelse(cor.m[id1, id2]>0, "+", "-")
  tab <- cbind(rbind(tab1,tab4), rbind(tab3, tab2))
  diag(tab) <- n2tab(diag(ssQ[id,id]))
  if(is.null(nome)){
    colnames(tab) <- colnames(cor.m)
    rownames(tab) <- rownames(cor.m)
  }
  
  # drifts estimados
  id <- grepl("slope", colnames(x$out$alphahat))
  mu <- utils::tail(x$out$alphahat[,id], 1)
  
  return(list(ar = n2tab(arpars), mu = n2tab(mu), tab = tab, model=x$out$model))
}



#' Intervalo de confianca para ucmodel
#' 
#' Calcula o intervalo de confianca de Bonferroni para um
#' componente estimado por ucmodel
#' 
#' @param ucm um modelo estimado por ucmodel
#' @param state o nome do estado estimado para o qual se deve
#' construir o intervalo de confianca
#' @param type tipo do estado usado, suavizado (padrao) ou filtrado 
#' @param ic o nivel de confianca do intervalo, entre 0 e 1
#' 
#' @return uma lista com 
#' \item{\code{upper}}{limite superior do intervalo}
#' \item{\code{lower}}{limite inferior do intervalo}
#' 
#' @export
#'    

ic.ucmodel <- function(ucm, state, type = "smooth", ic = 0.95){
  if(ic >= 1 | ic <=0){
    stop("ic fora do intervalo 0 < ic < 1")
  }
  if(type == "smooth"){
    n <- dim(ucm$out$alphahat)[1]
    id <- grep(state, colnames(ucm$out$alphahat))
    if(length(id)==1){
      v <- apply(ucm$out$V, 3, function(x) x[id,id])
    }else{
      v <- apply(ucm$out$V, 3, function(x) diag(x[id,id]))
    }
    v <- matrix(v, ncol = n)
    tn <- abs(qt(((1-ic)/2), df = n))
    suppressWarnings(v <- apply(v, 2, sqrt))
    v[is.nan(v)] <- NA
    v <- matrix(v, ncol = n)
    upper <- ucm$out$alphahat[,id] + t(tn*v)
    lower <- ucm$out$alphahat[,id] - t(tn*v)
  }
  if(type == "filter"){
    n <- dim(ucm$out$a)[1]
    id <- grep(state, colnames(ucm$out$a))
    if(length(id)==1){
      v <- apply(ucm$out$P, 3, function(x) x[id,id])
    }else{
      v <- apply(ucm$out$P, 3, function(x) diag(x[id,id]))
    }
    v <- matrix(v, ncol = n)
    tn <- abs(qt(((1-ic)/2), df = n))
    suppressWarnings(v <- apply(v, 2, sqrt))
    v[is.nan(v)] <- NA
    v <- matrix(v, ncol = n)
    upper <- ucm$out$a[,id] + t(tn*v)
    lower <- ucm$out$a[,id] - t(tn*v)
    upper <- window(lag(upper, 1), start = start(ucm$out$alphahat), end = end(ucm$out$alphahat))
    lower <- window(lag(lower, 1), start = start(ucm$out$alphahat), end = end(ucm$out$alphahat))
  }
  v <- t(v)
  colnames(v) <- colnames(ucm$out$alphahat)[id]
  return(list(upper = upper, lower = lower, v = v))
}


#' Tabela com as estimativas do ucmodel
#' 
#' Tabela com as estimativas do ucmodel 
#' 
#' @param ucm objeto ucmodel com modelo estimado por ucmodel
#' @return uma tabela com as estimativas
#' 
#' @export
#' 

tab.ucmodel <- function(ucm){
  nomes <- colnames(ucm$fit$model$y)
  n <- length(nomes)
  pars <- pars.ucmodel(ucm)
  tab <- data.frame(parametro = c(paste("$\\sigma_{\\eta}_{", 1:n ,"}$"),
                                  paste("$\\sigma_{\\varepsilon_{", 1:n ,"}}$"),
                                  paste("$\\mu_{", 1:n ,"}$"),
                                  as.vector(apply(matrix(paste("_{", 1:n,"}$")), 1, function(x) paste("$\\phi_{", 1:n ,"}", x, sep=""))),
                                  "log verossimilhanca",
                                  "Correlacoes",
                                  as.vector(apply(matrix(paste("$\\rho_{\\eta_{", 1:(n-1), "}")), 1, function(x) paste(x, "\\eta_{", 2:n ,"}}$"))),
                                  as.vector(apply(matrix(paste("$\\rho_{\\eta_{", 1:n, "}")), 1, function(x) paste(x, "\\varepsilon_{", 1:n ,"}}$"))),
                                  as.vector(apply(matrix(paste("$\\rho_{\\varepsilon_{", 1:(n-1), "}")), 1, function(x) paste(x, "\\varepsilon_{", 2:n ,"}}$")))
                                  ),
                    valor = c(diag(pars$tab),
                              pars$mu,
                              as.vector(pars$ar),
                              n2tab(ucm$out$logLik),
                              " ",
                              t(pars$tab)[lower.tri(pars$tab)])
                    )
  return(tab)
  
}



#' Plota ucmodel com ggplot2
#' 
#' @param ucm um modelo estimado por ucmodel
#' @param state o nome do estado estimado para o qual se deve
#' plotar o grafico
#' @param type tipo do estado usado, suavizado (padrao) ou filtrado
#' @param ic o nivel de confianca do intervalo, entre 0 e 1
#' 
#' @return ggplot das series
#' 
#' @import ggplot2 zoo
#' 
#' @export
#' 

ucplot <- function(ucm, state, type = "smooth", ic = 0.95){
  if(type == "smooth"){
    id <- grep(state, colnames(ucm$out$alphahat))
    x <- zoo::as.zoo(ucm$out$alphahat[,id])
    y <- zoo::as.zoo(ucm$out$model$y)
  }
  if(type == "filter"){
    id <- grep(state, colnames(ucm$out$a))
    x <- zoo::as.zoo(ucm$out$a[,id])
    y <- zoo::as.zoo(ucm$out$model$y)
    x <- window(lag(x, 1), start = start(y), end = end(y))
  }

  
  df.x <- zoo::fortify.zoo(x, melt = TRUE)
  df.y <- zoo::fortify.zoo(y, melt = TRUE)
  ic <- ic.ucmodel(ucm, state = state, type = type, ic = ic)
  ic.upper <- zoo::fortify.zoo(ic$upper, melt = TRUE)
  ic.lower <- zoo::fortify.zoo(ic$lower, melt = TRUE)
  df <- ggplot2::fortify(cbind(df.x, Value2=df.y[,3], upper=ic.upper[,3], lower=ic.lower[,3]), index.name = "Index")
  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = Index, y = Value))
  p <- p + ggplot2::geom_line(data = df, ggplot2::aes(x = Index, y = Value2),
                              linetype="dotted", size = 1/2, colour="red")
  p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper),
                                alpha=0.3)
  p <- p + ggplot2::geom_line(size = 1/2, alpha = 3/4)
  p <- p + ggplot2::facet_grid(Series ~ ., scales = "free_y") 
  #p <- p + ggplot2::facet_wrap(~ Series, scales = "free_y")
  p <- p + ggplot2::labs(y="", x="")
  p <- p + ggplot2::theme_bw()
  return(p)  
}
