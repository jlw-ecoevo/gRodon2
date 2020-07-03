
#' Get Codon Statistics
#'
#' Wrapper to estimate codon statistics  for an organism (used to fit gRodon).
#' It is assumed that gene names contain annotations of ribosomal proteins.
#'
#' @param gene_file  path to CDS-containing fasta file
getStatistics <- function(gene_file){
  genes <- readDNAStringSet(gene_file)
  highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)
  codon_stats <- getCodonStatistics(genes, highly_expressed)
  return(as.list(codon_stats))
}

#' Get Codon Statistics for All Genomes In a Directory
#'
#' This function gets the codon usage statistics for all CDS files in a directory (used to fit gRodon).
#' It is assumed that gene names contain annotations of ribosomal proteins.
#'
#' @param directory path to directory containing annotated CDS files
getStatisticsBatch <- function(directory, mc.cores = 1){
  gene_files <- list.files(directory)
  gene_paths <- paste0(directory,gene_files)
  cu <- mclapply(X = gene_paths, FUN = getStatistics, mc.cores = mc.cores) %>%
    do.call("rbind", .) %>%
    as.data.frame(stringsAsFactors = FALSE) %>%
    mutate(File=gene_files)
  return(cu)
}

#' Fit gRodon models
#'
#' This function fits the gRodon models
#'
#' @param stat_data dataframe with codon usage statistics and known doubling times
fitModels <- function(stat_data, stat_data_extremo){
  bc_milc <- boxcox(d~CUBHE+ConsistencyHE+CPB,data=stat_data)
  lambda_milc <- bc_milc$x[which.max(bc_milc$y)]

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

  return(list(gRodon_model_base,
              gRodon_model_temp,
              gRodon_model_partial,
              gRodon_model_partial_temp,
              gRodon_model_meta,
              gRodon_model_meta_temp,
              lambda_milc))
}
