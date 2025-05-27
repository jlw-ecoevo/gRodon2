## JLW 2020 - Fit gRodon

# Load Packages ----------------------------------------------------------------

library(dplyr)
library(gRodon)
library(parallel)
library(Biostrings)
library(coRdon)
library(matrixStats)


# Calculate Codon Usage Statistics ---------------------------------------------

genomes_path <- ""
trimlens <- seq(120,300,10)

cu <- list()
cu_i <- list()
for(i in 1:length(trimlens)){
  print(paste("truncate forward: ",trimlens[i]))
  cu[[i]] <- gRodon:::getStatisticsBatch(genomes_path,
                                         mc.cores = 28,
                                         fragments = TRUE,
                                         bg = "all",
                                         trimlen = trimlens[i],
                                         trimside = "start")
  cu_i[[i]] <- gRodon:::getStatisticsBatch(genomes_path,
                                         mc.cores = 28,
                                         fragments = TRUE,
                                         bg = "individual",
                                         trimlen = trimlens[i],
                                         trimside = "start")
}

cu_start <- list(cu=cu,cu_i=cu_i,trimlens=trimlens)

setwd("~/gRodon/inst/extdata/")
save(cu_start, file = "CodonStatistics_truncated_genes_start.rda")

genomes_path <- ""
trimlens <- seq(120,300,10)

cu <- list()
cu_i <- list()
for(i in 1:length(trimlens)){
  print(paste("truncate backward: ",trimlens[i]))
  cu[[i]] <- gRodon:::getStatisticsBatch(genomes_path,
                                         mc.cores = 28,
                                         fragments = TRUE,
                                         bg = "all",
                                         trimlen = trimlens[i],
                                         trimside = "end")
  cu_i[[i]] <- gRodon:::getStatisticsBatch(genomes_path,
                                           mc.cores = 28,
                                           fragments = TRUE,
                                           bg = "individual",
                                           trimlen = trimlens[i],
                                           trimside = "end")
}

cu_end <- list(cu=cu,cu_i=cu_i,trimlens=trimlens)

setwd("~/gRodon/inst/extdata/")
save(cu_end, file = "CodonStatistics_truncated_genes_end.rda")

