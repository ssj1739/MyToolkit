# Neural net for protein prediction
library(nnet)
library(dplyr)
library(tidyr)
library(magrittr)

load("~/Downloads/RNAseq.RData")
load("~/Downloads/Proteins20180119.rdata")

df = gather(RNAseq$dat %>% set_colnames(RNAseq$gene))

match(df$key, RNAseq$ge)