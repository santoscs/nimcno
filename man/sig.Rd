% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/tabelas.R
\name{sig}
\alias{sig}
\title{Significancia de um teste}
\usage{
sig(stat, cval, comp = TRUE)
}
\arguments{
\item{stat}{valor da estatistica do teste}

\item{cval}{vetor com os valores criticos do teste para 10\%, 5\% e 1\%}

\item{comp}{logical. If TRUE (the default) the statistics is print too}
}
\value{
um objeto "character" com valor da estatistica do teste formatado
junto com "***" para 1\%, "**" para 5\%, "*" para 10\% e " " nao significativo
}
\description{
Fornece o valor da estatistica de um teste junto com a indicacao
da significancia para valores 10\%, 5\% e 1\%.
}
