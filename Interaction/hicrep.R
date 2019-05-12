#Description: Code for run hicrep to calculate the stratum-adjusted correlation coefficient (scc). Please see https://bioconductor.org/packages/release/bioc/manuals/hicrep/man/hicrep.pdf for details


#/usr/bin/Rscript
library("hicrep")
options(echo=TRUE)
args<-commandArgs(trailingOnly=TRUE)
input1<-paste(args[1])
input2<-paste(args[2])
HiCR1<-read.table(args[1])
HiCR2<-read.table(args[2])
                                                            
#Estimate the optimial smoothing neighborhood size parameter
h_hat <- htrain(HiCR1, HiCR2, 1000000, 5000000, 0:10)
h_hat <- 0
processed <- prep(HiCR1, HiCR2, 1000000, 1, 5000000)
scc.out <- get.scc(processed, 1000000, 5000000)
scc.out$scc
scc.out$std
