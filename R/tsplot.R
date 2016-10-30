#' Plota series temporais com ggplot2
#' 
#' @param x objeto ts ou mts com as series temporais
#' @param y (opicional) objeto ts ou mts com dimensao de x para ser
#' plotado junto com x no mesmo grafico 
#' @param escala Are scales shared across all facets
#'  ("fixed"), or do they vary across 
#'  rows ("free_x"), columns (the default, "free_y"), or both 
#'  rows and columns ("free")
#' @param facet as series sao plotadas em graficos diferente (facet = TRUE, the default),
#' ou no mesmo grafico (facet = FALSE)
#' @param name optional name for ts univariate
#' 
#' @return ggplot das series
#' 
#' @import ggplot2 zoo
#' 
#' @export
#' 

tsplot <- function(x, y = NULL, escala = 'free_y', facet = TRUE, name = NULL){
  nseries <- NCOL(x)
  ntime <- NROW(x)
  x <- zoo::as.zoo(x)
  df.x <- zoo::fortify.zoo(x, melt = TRUE)
  if(nseries==1 & !is.null(name)){
    df.x[,"Series"] <- rep(name, ntime)
  }
  if(!is.null(y)){
    y <- zoo::as.zoo(y)
    df.y <- zoo::fortify.zoo(y, melt = TRUE)
    if(facet){
      df <- ggplot2::fortify(cbind(df.x, Value2=df.y[,3]), index.name = "Index")
      p <- ggplot2::ggplot(data = df, ggplot2::aes(x = Index, y = Value))
      p <- p + ggplot2::geom_line(data = df, ggplot2::aes(x = Index, y = Value2),
                                  linetype="dotted", size = 1/2)
      p <- p + ggplot2::geom_line(size = 1/2, alpha = 3/4)  
      p <- p + ggplot2::facet_grid(Series ~ ., scales = "free_y") 
      #p <- p + ggplot2::facet_wrap(~ Series, scales = "free_y")
    }else{
      p <- ggplot(df.x, aes(x = Index, y = Value))
      p <- p + geom_line(data = df.y, aes(x = Index, y = Value, group = Series), size = 1/2, alpha = 3/4)
      p <- p + geom_line(linetype="dotted", size = 1/2)  # Drawing the "overlayer"
    }
    p <- p + ggplot2::labs(y="", x="")
    p <- p + ggplot2::theme_bw()
    return(p)  
  }
  if(!facet){
    p <- ggplot2::ggplot(data = df.x, ggplot2::aes(x = Index, y = Value, color=Series))
    p <- p + ggplot2::geom_line(size = 1/2, alpha = 3/4)  
    p <- p + ggplot2::theme_bw()
    return(p)  
    
  }else{
    p <-ggplot2::ggplot(df.x, ggplot2::aes(x=Index, y=Value, group_by())) +
      ggplot2::geom_line(size = 1/2, alpha = 3/4) +
      ggplot2::facet_grid(Series ~ ., scales = escala) +
      ggplot2::labs(y="", x="") +
      ggplot2::theme_bw()
  }
  return(p)
}


#' Intervalo de confianca 
#' 
#' Calcula o intervalo de confianca para uma estimativa pontual
#' 
#' @param x estimativa pontual
#' @param dp desvio padrao da estimativa pontual
#' @param nc nivel de confianca do intervalo, entre 0 e 1
#' 
#' @return uma lista com 
#' \item{\code{upper}}{limite superior do intervalo}
#' \item{\code{lower}}{limite inferior do intervalo}
#' 
#' @export
#'    

ic <- function(x, dp, nc = 0.95){
  if(nc >= 1 | nc <=0){
    stop("nc fora do intervalo 0 < nc < 1")
  }
  n <- dim(x)[1]
  # margem de erro
  tn <- abs(qt(((1-nc)/2), df = n))
  upper <- x + tn*dp
  lower <- x - tn*dp
  return(list(upper = upper, lower = lower))
}

#' Plot com ggplot2
#' 
#' @param x series a serem plotadas
#' @param y opcional serie ou series a serem plotadas juntas para comparacao
#' @param upper limite superior para regiao sombreada
#' @param lower limite inferior para regiao sombreada
#' 
#' @return ggplot das series
#' 
#' @import ggplot2 zoo
#' 
#' @export
#' 

biplot <- function(x, y, upper, lower){
  df.x <- zoo::fortify.zoo(x, melt = TRUE)
  df.y <- zoo::fortify.zoo(y, melt = TRUE)
  ic.upper <- zoo::fortify.zoo(upper, melt = TRUE)
  ic.lower <- zoo::fortify.zoo(lower, melt = TRUE)
  df <- ggplot2::fortify(cbind(df.x, Value2=df.y[,3], upper=ic.upper[,3], lower=ic.lower[,3]), index.name = "Index")
  p <- ggplot2::ggplot(data = df, ggplot2::aes(x = Index, y = Value))
  p <- p + ggplot2::geom_line(data = df, ggplot2::aes(x = Index, y = Value2),
                              linetype="dotted", size = 1/2)
  p <- p + ggplot2::geom_ribbon(ggplot2::aes(ymin=lower, ymax=upper),
                                alpha=0.3)
  p <- p + ggplot2::geom_line(size = 1/2, alpha = 3/4)
  p <- p + ggplot2::facet_grid(Series ~ ., scales = "free_y") 
  p <- p + ggplot2::labs(y="", x="")
  p <- p + ggplot2::theme_bw()
  return(p)  
}



