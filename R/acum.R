#' Acumula series em percentual ao mes em m meses
#' 
#' Transforma uma serie mensal dada em percentual ao mes 
#' em uma serie mensal com percentual nos ultimos m meses
#'
#' @param x A time series univariate
#' @param m number of monthes
#' 
#' @return A time series univariate 
#' 
#' @import zoo stats
#' 
#' @export

acum<-function(x, m=12){
  # input:
  # x(ts): serie a ser acumulada
  #output: 
  # x12(ts): serie acumulada
  
  x <- zoo::as.zoo(x)
  x12 <- zoo::rollapplyr(x, width=m, function(x) (prod(1+x/100)-1)*100)
  x12 <- stats::as.ts(x12)
  return(x12)
}
