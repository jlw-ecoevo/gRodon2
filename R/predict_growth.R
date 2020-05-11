
#' Predict Growth Rate
#'
#' This function predicts the growth rate of a prokaryotic organism based on
#' genomic codon usage patterns (Weissman et. al. TBD).
#'
#' @param genes DNAStringSet object holding all nucleotide sequences from a
#' genome.
#' @param highly_expressed Logical vector describing the set of highly expressed
#' genes. Must be of same lend as \code{genes}.
#' @param mode Whether to run prediction in full, partial, or metagenome mode
#' (by default gRodon applies the full model)
#' @param temperature Optimal growth temperature. By default this is set as
#' "none" and we do not guarantee good results for non-mesophilic organisms.
#' @return gRodon returns a list with the following elements:
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
#' # Load in example genome (Streptococcus pyogenes M1, downloaded from RefSeq)
#' # included with gRodon
#' path_to_genome <- system.file('extdata',
#'   'GCF_000349925.2_ASM34992v2_cds_from_genomic.fna.gz',
#'   package = 'gRodon')
#' genes <- readDNAStringSet(path_to_genome)
#'
#' # Search pre-existing annotations for ribosomal proteins, which we
#' # will use as our set of highly expressed genes
#' highly_expressed <- grepl("ribosomal protein",names(genes),ignore.case = T)
#'
#' # Run the gRodon growth prediction pipeline
#' predictGrowth(genes, highly_expressed)
#'
#' # Run gRodon with temperature option (not needed for mesophiles, gRodon not
#' # validated on extremophiles, use with care)
#' predictGrowth(genes, highly_expressed, temperature = 37)
#'
#' @export
#' @importFrom magrittr "%>%"
#' @importFrom dplyr group_by mutate count summarise
#' @import coRdon
#' @import Biostrings
predictGrowth <- function(genes,
                          highly_expressed,
                          mode = "full",
                          temperature = "none"){

  if(sum(highly_expressed)<10){
    warning("Less than 10 highly expressed genes provided, performance may suffer")
  }

  # Calculate codon data
  codon_stats <- getCodonStatistics(genes, highly_expressed)

  # Predict growth rate (stored models - sysdata.rda)
  if(temperature == "none" & mode=="full"){
    pred <- stats::predict.lm(gRodon_model_base,
                              newdata = codon_stats,
                              interval = "confidence")

  } else if(temperature == "none" & mode=="metagenome"){
    pred <- stats::predict.lm(gRodon_model_meta,
                              newdata = codon_stats,
                              interval = "confidence")

  } else if(temperature == "none" & mode=="partial"){
    pred <- stats::predict.lm(gRodon_model_partial,
                              newdata = codon_stats,
                              interval = "confidence")

  } else if(temperature != "none" & mode=="full"){
    codon_stats$OGT <- temperature
    pred <- stats::predict.lm(gRodon_model_temp,
                              newdata = codon_stats,
                              interval = "confidence")

  } else if(temperature != "none" & mode=="metagenome"){
    codon_stats$OGT <- temperature
    pred <- stats::predict.lm(gRodon_model_meta_temp,
                              newdata = codon_stats,
                              interval = "confidence")

  } else if(temperature != "none" & mode=="partial"){
    codon_stats$OGT <- temperature
    pred <- stats::predict.lm(gRodon_model_partial_temp,
                              newdata = codon_stats,
                              interval = "confidence")
  }

  #Transform back from box-cox
  pred_back_transformed <- boxcoxTransform(pred,
                                           lambda_milc,
                                           back_transform = TRUE)
  #attach prediction
  codon_stats$d <- pred_back_transformed[,"fit"]
  codon_stats$LowerCI <- pred_back_transformed[,"lwr"]
  codon_stats$UpperCI <- pred_back_transformed[,"upr"]

  #Return
  if(is.na(pred_back_transformed[,"fit"]) & pred[,"fit"]>6){
    warning("Estimated doubling time very long. Essentially goes to infinity
(gives NA value after back-transforming from box-cox).")
  } else if(pred_back_transformed[,"fit"]>5){
    warning("Estimated doubling time >5 hours. CUB signal saturates at approx. 5 hrs ...
gRodon may underestimate doubling times above this range.
Consider simply reporting as '>5hrs'.
(In other words, this microbe definitely grows slowly, but we can't tell you quite how slowly).")
  }
  return(as.list(codon_stats))
}
