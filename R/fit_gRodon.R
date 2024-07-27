
#' Get Codon Statistics
#'
#' Wrapper to estimate codon statistics  for an organism (used to fit gRodon).
#' It is assumed that gene names contain annotations of ribosomal proteins.
#'
#' @param gene_file  path to CDS-containing fasta file
getStatistics <- function(gene_file, genetic_code = "11", bg = "all"){
  print(gene_file)
  genes <- readDNAStringSet(gene_file)
  highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)
  # print(table(highly_expressed))
  if(sum(highly_expressed)<10){
    return(NULL)
  } else {
    if(bg=="all"){
      codon_stats <- try(getCodonStatistics(genes,
                                            highly_expressed,
                                            genetic_code = genetic_code))
    } else if(bg=="individual"){
      codon_stats <- try(getCodonStatistics_i(genes,
                                            highly_expressed,
                                            genetic_code = genetic_code))
    } else{
      stop("Feature in testing, please set bg==\"all\" for normal gRodon behavior")
    }
    codon_stats[["File"]] <- basename(gene_file)
    if(!inherits(codon_stats,"try-error")){
      return(as.list(codon_stats))
    } else {
      return(NULL)
    }
  }
}

#' Get Codon Statistics for All Genomes In a Directory
#'
#' This function gets the codon usage statistics for all CDS files in a directory (used to fit gRodon).
#' It is assumed that gene names contain annotations of ribosomal proteins.
#'
#' @param directory path to directory containing annotated CDS files
getStatisticsBatch <- function(directory, genetic_code="11", mc.cores = 1, bg = "all"){
  gene_files <- list.files(directory)
  gene_paths <- paste0(directory,gene_files)
  cu <- mclapply(X = gene_paths,
                 FUN = getStatistics,
                 genetic_code = genetic_code,
                 bg = bg,
                 mc.cores = mc.cores) %>%
    do.call("rbind", .) %>%
    as.data.frame(stringsAsFactors = FALSE) #%>%
    # dplyr::mutate(File=gene_files)
  return(cu)
}


getStatisticsBatchWindows <- function(directory, genetic_code="11", mc.cores = 1, bg = "all"){
  gene_files <- list.files(directory)
  gene_paths <- paste0(directory,gene_files)
  cl <- makeCluster(mc.cores)
  clusterEvalQ(cl, c(library(dplyr),
                     library(Biostrings),
                     library(coRdon),
                     library(gRodon),
                     library(matrixStats)))
  clusterExport(cl=cl, as.list(ls()),
                envir=environment())
  cu <- parLapply(cl = cl,
                  X = gene_paths,
                  fun = getStatistics,
                  genetic_code = genetic_code,
                  bg = bg)
  stopCluster(cl)
  # cu <- cu %>%
  #   do.call("rbind", .) %>%
  #   as.data.frame(stringsAsFactors = FALSE)
  return(cu)
}



#' Fit gRodon models
#'
#' This function fits the gRodon models
#'
#' @param stat_data dataframe with codon usage statistics and known doubling times
fitModels <- function(stat_data, stat_data_extremo){
  bc_milc <- boxcox(d~CUBHE+ConsistencyHE+CPB,data=stat_data)
  bc_milc_euk <- boxcox(d~CUBHE,data=stat_data)
  lambda_milc <- bc_milc$x[which.max(bc_milc$y)]
  lambda_milc_euk <- bc_milc_euk$x[which.max(bc_milc_euk$y)]

  #Full gRodon
  gRodon_model_base <-
    lm(boxcoxTransform(d, lambda_milc) ~ CUBHE+ConsistencyHE+CPB,data=stat_data)
  gRodon_model_temp <-
    lm(boxcoxTransform(d, lambda_milc) ~ CUBHE+ConsistencyHE+CPB+OGT,data=stat_data_extremo)

  # Partial genome mode
  gRodon_model_partial <-
    lm(boxcoxTransform(d, lambda_milc) ~ CUBHE+ConsistencyHE,data=stat_data)
  gRodon_model_partial_temp <-
    lm(boxcoxTransform(d, lambda_milc) ~ CUBHE+ConsistencyHE+OGT,data=stat_data_extremo)

  # Metagenome mode
  gRodon_model_meta <-
    lm(boxcoxTransform(d, lambda_milc) ~ CUBHE,data=stat_data)
  gRodon_model_meta_temp <-
    lm(boxcoxTransform(d, lambda_milc) ~ CUBHE+OGT,data=stat_data_extremo)

  # Eukaryotic mode
  gRodon_model_euk <-
    lm(boxcoxTransform(d, lambda_milc_euk) ~ CUBHE,data=stat_data)
  gRodon_model_euk_temp <-
    lm(boxcoxTransform(d, lambda_milc_euk) ~ CUBHE+OGT,data=stat_data_extremo)

  return(list(gRodon_model_base,
              gRodon_model_temp,
              gRodon_model_partial,
              gRodon_model_partial_temp,
              gRodon_model_meta,
              gRodon_model_meta_temp,
              lambda_milc,
              gRodon_model_euk,
              gRodon_model_euk_temp,
              lambda_milc_euk))
}


#' Fit gRodon GC-corrected models
#'
#' This function fits the gRodon GC-corrected metagenome mode models
#'
#' @param stat_data dataframe with codon usage statistics and known doubling times
fitGCModels <- function(stat_data, stat_data_extremo){

  stat_data <- stat_data %>%
    mutate(dCUB=(CUB-CUBHE)/CUB)
  stat_data_extremo <- stat_data_extremo %>%
    mutate(dCUB=(CUB-CUBHE)/CUB)

  bc_meta <- boxcox(d~dCUB+GCdiv,data=stat_data)
  lambda_meta <- bc_meta$x[which.max(bc_meta$y)]

  #new metagenome mode
  meta_model_base <-
    lm(boxcoxTransform(d, lambda_meta) ~ dCUB+GCdiv,data=stat_data)
  meta_model_temp <-
    lm(boxcoxTransform(d, lambda_meta) ~ dCUB+GCdiv+OGT,data=stat_data_extremo)

  #new metagenome mode without gc cor
  meta_nogc_model_base <-
    lm(boxcoxTransform(d, lambda_meta) ~ dCUB,data=stat_data)
  meta_nogc_model_temp <-
    lm(boxcoxTransform(d, lambda_meta) ~ dCUB+OGT,data=stat_data_extremo)

  #old metagenome mode (for individual-gene fitting)
  bc_oldmeta <- boxcox(d~CUBHE,data=stat_data)
  lambda_oldmeta <- bc_oldmeta$x[which.max(bc_oldmeta$y)]
  meta_old_model_base <-
    lm(boxcoxTransform(d, lambda_oldmeta) ~ CUBHE,data=stat_data)
  meta_old_model_temp <-
    lm(boxcoxTransform(d, lambda_oldmeta) ~ CUBHE+OGT,data=stat_data_extremo)

  return(list(meta_model_base,
              meta_model_temp,
              meta_nogc_model_base,
              meta_nogc_model_temp,
              lambda_meta,
              meta_old_model_base,
              meta_model_temp,
              lambda_oldmeta))
}
