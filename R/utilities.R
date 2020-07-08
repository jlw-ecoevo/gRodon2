

filterSeq <- function(genes, highly_expressed, length_threshold = 240, depth_of_coverage = NULL){
  # Gene lengths >= 80 codons for good behavior (see coRdon documentation)
  # Throw out any genes w/ ambiguous frame
  pass_filter <- (width(genes) > length_threshold) & ((width(genes) %% 3) == 0)
  highly_expressed <- highly_expressed[pass_filter]
  genes <- genes[pass_filter]
  if(!is.null(depth_of_coverage)){
    depth_of_coverage <- depth_of_coverage[pass_filter]
  }

  # Remove genes containing ambiguous bases
  ambiguous_bases <- grepl("K|M|R|S|W|Y|N|V|H|D|B",genes)
  highly_expressed <- highly_expressed[!ambiguous_bases]
  genes <- genes[!ambiguous_bases]
  if(!is.null(depth_of_coverage)){
    depth_of_coverage <- depth_of_coverage[!ambiguous_bases]
  }

  # Warn user if genes have been filtered
  if(sum(!pass_filter) > 0){
    warning(paste("There were", sum(!pass_filter),
                  "genes either with lengths not multiples of 3 or not above length threshold (default 240bp), these genes have been ignored"))
  }
  if(sum(ambiguous_bases) > 0){
    warning(paste("There were", sum(ambiguous_bases),
                  "genes with ambiguous bases, these genes have been ignored"))
  }
  if(sum(highly_expressed) == 0){
    stop("No highly expressed genes after filtering, unable to compute growth rate")
  }

  return(list(Genes=genes,
              HE = highly_expressed,
              Filtered = sum(ambiguous_bases)+sum(!pass_filter),
              Depth = depth_of_coverage))
}

boxcoxTransform <- function(x, lambda, back_transform = F) {
  if (back_transform == TRUE) {
    (x*lambda +1)^(1/lambda)  %>% return()
  } else {
    (((x^lambda) - 1) / lambda) %>% return()
  }
}
