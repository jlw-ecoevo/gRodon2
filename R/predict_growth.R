
#' Predict Growth Rate
#'
#' This function predicts the growth rate of a prokaryotic organism based on
#' genomic codon usage patterns (Weissman et. al. TBD).
#'
#' @param genes DNAStringSet object holding all nucleotide sequences from a
#' genome.
#' @param highly_expressed Logical vector describing the set of highly expressed
#' genes. Must be of same lend as \code{genes}.
#' @param metagenome Whether to run prediction in metagenome mode (no codon pair
#' bias)
#' @param temperature Optimal growth temperature. By default this is set as
#' "none" and we do not guarantee good results for non-mesophilic organisms.
#' @return A list with the following components is returned:
#' \describe{
#'   \item{CUBHE}{Median codon usage bias of the highly expressed genes (MILC)
#'   calculated using the genome-wide codon usage as the expected bias}
#'   \item{ConsistencyHE}{Median codon usage bias of the highly expressed genes
#'   (MILC) calculated using the codon usage of highly expressed genes as the
#'   expected bias}
#'   \item{CPB}{Genome-wide codon pair bias (Coleman et al. 2008)}
#'   \item{FilteredSequences}{Number of gene sequences filtered out during
#'   calulation (due to length and/or presence of ambiguous bases)}
#'   \item{d}{Predicted doubling time}
#'   \item{LowerCI}{Lower CI of \code{d} (2.5%)}
#'   \item{UpperCI}{Upper CI of \code{d} (97.5%)}
#' }
#' @examples
#' genes <- readDNAStringSet(path_to_genome)
#' highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)
#' getGrowth(genes, highly_expressed)
#' @export
#' @import dplyr
#' @import coRdon
#' @import Biostrings
predictGrowth <- function(genes,
                      highly_expressed,
                      metagenome = FALSE,
                      temperature = "none"){

  # Calculate codon data
  codon_stats <- getCodonStatistics(genes, highly_expressed)

  # Predict growth rate (stored models - sysdata.rda)
  if(temperature == "none" & metagenome == FALSE){
    pred <- stats::predict.lm(gRodon_model_base,
                              newdata = codon_stats,
                              interval = "confidence")

  } else if(temperature == "none" & metagenome == TRUE){
    pred <- stats::predict.lm(gRodon_model_meta,
                              newdata = codon_stats,
                              interval = "confidence")

  } else if(temperature != "none" & metagenome == FALSE){
    pred <- stats::predict.lm(gRodon_model_temp,
                              newdata = codon_stats,
                              interval = "confidence")

  } else if(temperature != "none" & metagenome == TRUE){
    pred <- stats::predict.lm(gRodon_model_meta_temp,
                              newdata = codon_stats,
                              interval = "confidence")
  }

  #Transform back from box-cox
  pred_back_transformed <- boxcoxTransform(pred,
                                           lambda_milc,
                                           back_transform = TRUE)

  #return prediction
  codon_stats$d <- pred_back_transformed[,"fit"]
  codon_stats$LowerCI <- pred_back_transformed[,"lwr"]
  codon_stats$UpperCI <- pred_back_transformed[,"upr"]
  return(as.list(codon_stats))
}
