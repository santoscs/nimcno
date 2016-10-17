### 2 transforma os dados  ###
### Projeto: nimcno        ###
##############################

## importa os dados
rm(list = ls()) # limpa workspace

dados <- read.csv2("data-raw/dados_95_16.csv", na.strings="-")
# log m1
dados <- dplyr::mutate(dados, logm1=log(m1))
# pib real
dados <- dplyr::mutate(dados, pibr= pib-(igp.m/100)*pib)
# log do pib real
dados <- dplyr::mutate(dados, logpibr=log(pibr))
# selic real ex-post(igp-m)
dados <- dplyr::mutate(dados, selicr=(((1 + selic/100)/(1+igp.m/100))-1)*100)
# exclui as series não usadas

tsm <- ts(dados[,-1], start = c(1995,1), frequency = 12)

# # ajuste sazonal
# s <- function(x){
#   #install.packages(c("seasonal"))
#   seasonal::seas(x)$data[,"final"]
# }
# 
# #log m1 ajuste sazonal
# tsm[,"logm1"] <- s(tsm[,"logm1"])
# 
# #log m1 ajuste sazonal
# tsm[,"logpibr"] <- s(tsm[,"logpibr"])




## Salva os dados no pacote
macro95 <- tsm
devtools::use_data(macro95, overwrite = TRUE)


write.csv2(macro95, file = "data-raw/dados95-16.csv")
