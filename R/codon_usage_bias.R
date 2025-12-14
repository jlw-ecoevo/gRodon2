
shuffleGenes <- function(gene){
  DNAString(stringi::stri_rand_shuffle(as.character(gene))) %>%
    return()
}

singleMILC <- function(gene,n=100,method="MILC"){
  bg <- replicate(n,shuffleGenes(gene)) %>% DNAStringSet()
  # print(bg)
  # print(gene)
  genes_bg <- c(DNAStringSet(gene),bg)
  l1 <- logical(n+1)
  l1[1] <- T

  if(method=="MILC"){
    MILC(codonTable(genes_bg))[1,1] %>%
      return()

  } else if(method == "ENCprime"){
    ENCprime(codonTable(genes_bg))[1,1] %>%
      return()

  } else if(method == "B"){
    B(codonTable(genes_bg))[1,1] %>%
      return()

  } else if(method == "SCUO"){
    SCUO(codonTable(genes_bg))[1] %>%
      return()

  } else if(method == "MCB"){
    MCB(codonTable(genes_bg))[1,1] %>%
      return()

  } else {
    stop("Error: Please pick an implemented method (MILC, MILCgenomic, consistency)")
  }

}


CUBi <- function(genes,highly_expressed,n_le=100,method="MILC"){
  genes_list <- as.list(as.character(genes)) %>%
    lapply(DNAString)
  genes_list_HE <- genes_list[highly_expressed]
  genes_list_LE <- genes_list[!highly_expressed][sample(1:sum(!highly_expressed),n_le,replace=T)]

  x_HE <- lapply(genes_list_HE,singleMILC,method=method) %>% unlist()
  x_LE <- lapply(genes_list_LE,singleMILC,method=method) %>% unlist()

  return(c(median(x_LE),median(x_HE)))
}

getCUB <- function(fna_tab, highly_expressed, method = "MILC", genetic_code = "11"){

  if(method == "MILC"){
    # MILC estimate of codon usage bias from coRdon
    #   - Bias of highly expessed genes, with genome-wide codon usage considered
    #     expected codon usage
    x <- MILC(fna_tab,id_or_name2 = genetic_code)
    # print(x[highly_expressed, 1])
    x[highly_expressed, 1] %>% median() %>% return()

  } else if(method == "MILCgenomic"){
    # MILC estimate of codon usage bias from coRdon
    #   - Bias of all genes, with genome-wide codon usage considered
    #     expected codon usage
    x <- MILC(fna_tab,id_or_name2 = genetic_code)
    x[, 1] %>% median() %>% return()

  } else if(method == "consistency"){
    # Consistency of CUB across highly expressed genes
    #   - Bias of highly expessed genes, with codon usage of highly expressed
    #     gene set considered expected codon usage
    milc <- MILC(fna_tab,
                 subsets = list(HE = highly_expressed),
                 id_or_name2 = genetic_code)
    milc[highly_expressed, 2] %>% mean() %>% return()

  } else if(method == "ENCprime"){
      # ENCprime estimate of codon usage bias from coRdon
      #   - Bias of highly expessed genes, with genome-wide codon usage considered
      #     expected codon usage
      x <- ENCprime(fna_tab,id_or_name2 = genetic_code)
      # print(x[highly_expressed, 1])
      x[highly_expressed, 1] %>% median() %>% return()

  } else if(method == "B"){
    # B estimate of codon usage bias from coRdon
    #   - Bias of highly expessed genes, with genome-wide codon usage considered
    #     expected codon usage
    x <- B(fna_tab,id_or_name2 = genetic_code)
    # print(x[highly_expressed, 1])
    x[highly_expressed, 1] %>% median() %>% return()

  } else if(method == "SCUO"){
    # SCUO estimate of codon usage bias from coRdon
    #   - Bias of highly expessed genes, with genome-wide codon usage considered
    #     expected codon usage
    x <- SCUO(fna_tab,id_or_name2 = genetic_code)
    # print(x[highly_expressed, 1])
    x[highly_expressed] %>% median() %>% return()

  } else if(method == "MCB"){
    # MCB estimate of codon usage bias from coRdon
    #   - Bias of highly expessed genes, with genome-wide codon usage considered
    #     expected codon usage
    x <- MCB(fna_tab,id_or_name2 = genetic_code)
    # print(x[highly_expressed, 1])
    x[highly_expressed, 1] %>% median() %>% return()

  } else {
    stop("Error: Please pick an implemented method (MILC, MILCgenomic, consistency)")
  }
}

getWeightedCUB <- function(fna_tab, highly_expressed, depth_of_coverage, genetic_code = "11"){
    x <- MILC(fna_tab,id_or_name2 = genetic_code)
    x[, 1] %>%
      weightedMedian(., w = depth_of_coverage) %>%
      return()
}

getWeightedCUBHE <- function(fna_tab, highly_expressed, depth_of_coverage, genetic_code = "11",method="MILC"){

  if(method=="MILC"){
    x <- MILC(fna_tab,id_or_name2 = genetic_code)
    x[highly_expressed, 1] %>%
      weightedMedian(., w = depth_of_coverage[highly_expressed]) %>%
      return()
  }

  else if(method == "ENCprime"){
    # MILC estimate of codon usage bias from coRdon
    #   - Bias of highly expessed genes, with genome-wide codon usage considered
    #     expected codon usage
    x <- ENCprime(fna_tab,id_or_name2 = genetic_code)
    # print(x[highly_expressed, 1])
    x[highly_expressed, 1] %>%
      weightedMedian(., w = depth_of_coverage[highly_expressed]) %>%
      return()

  } else if(method == "B"){
    # MILC estimate of codon usage bias from coRdon
    #   - Bias of highly expessed genes, with genome-wide codon usage considered
    #     expected codon usage
    x <- B(fna_tab,id_or_name2 = genetic_code)
    # print(x[highly_expressed, 1])
    x[highly_expressed, 1] %>%
      weightedMedian(., w = depth_of_coverage[highly_expressed]) %>%
      return()

  } else if(method == "SCUO"){
    # MILC estimate of codon usage bias from coRdon
    #   - Bias of highly expessed genes, with genome-wide codon usage considered
    #     expected codon usage
    x <- SCUO(fna_tab,id_or_name2 = genetic_code)
    # print(x[highly_expressed, 1])
    x[highly_expressed] %>%
      weightedMedian(., w = depth_of_coverage[highly_expressed]) %>%
      return()

  } else if(method == "MCB"){
    # MILC estimate of codon usage bias from coRdon
    #   - Bias of highly expessed genes, with genome-wide codon usage considered
    #     expected codon usage
    x <- MCB(fna_tab,id_or_name2 = genetic_code)
    # print(x[highly_expressed, 1])
    x[highly_expressed, 1] %>%
      weightedMedian(., w = depth_of_coverage[highly_expressed]) %>%
      return()

  } else {
    stop("Error: Please pick an implemented method (MILC)")
  }

}

getWeightedCUBi <- function(genes,highly_expressed,depth_of_coverage,n_le=100,method="MILC"){
  genes_list <- as.list(as.character(genes)) %>%
    lapply(DNAString)
  samp_le <- sample(1:sum(!highly_expressed),n_le,replace=T)
  genes_list_HE <- genes_list[highly_expressed]
  genes_list_LE <- genes_list[!highly_expressed][samp_le]

  x_HE <- lapply(genes_list_HE,singleMILC,method=method) %>% unlist()
  x_LE <- lapply(genes_list_LE,singleMILC,method=method) %>% unlist()

  wHE <- x_HE %>%
    weightedMedian(., w = depth_of_coverage[highly_expressed])
  wLE <- x_LE %>%
    weightedMedian(., w = depth_of_coverage[!highly_expressed][samp_le])

  return(c(median(wLE),median(wHE)))
}

getWeightedConsistency <- function(fna_tab, highly_expressed, depth_of_coverage, genetic_code = "11"){
  x <- MILC(fna_tab,
            subsets = list(HE = highly_expressed),
            id_or_name2 = genetic_code)
  sum(x[highly_expressed, 2]*(depth_of_coverage[highly_expressed]/sum(depth_of_coverage[highly_expressed]))) %>%
    return()
}

getCodonStatistics <- function(genes,
                               highly_expressed,
                               fragments = FALSE,
                               depth_of_coverage = NULL,
                               genetic_code = "11",
                               trimlen = NA,
                               trimside = "start",
                               all_metrics=F){

  if(sum(highly_expressed) == 0){
    stop("No highly expressed genes?")
  }

  #Remove short sequences and sequences that are not multiples of 3
  if(!is.na(fragments)){
    if(fragments == TRUE){
      filtered <- filterSeq(genes = genes,
                            highly_expressed = highly_expressed,
                            length_threshold = 120,
                            depth_of_coverage = depth_of_coverage)
    } else {
      filtered <- filterSeq(genes = genes,
                            highly_expressed = highly_expressed,
                            depth_of_coverage = depth_of_coverage)
    }
  } else {
    filtered <- filterSeq(genes = genes,
                          highly_expressed = highly_expressed,
                          length_threshold = 1,
                          depth_of_coverage = depth_of_coverage)
  }


  if(is.na(trimlen)){
    genes <- filtered$Genes
    highly_expressed <- filtered$HE
    depth_of_coverage <- filtered$Depth
  } else if (trimside == "start"){
    genes <- c(subseq(filtered$Genes[width(filtered$Genes)>=trimlen], 1, trimlen),
               filtered$Genes[width(filtered$Genes)<trimlen])
    ind_genes <- c(which(width(filtered$Genes)>=trimlen), which(width(filtered$Genes)<trimlen))
    highly_expressed <- filtered$HE[ind_genes]
    depth_of_coverage <- filtered$Depth[ind_genes]
  } else if (trimside == "end"){
    genes <- c(subseq(filtered$Genes[width(filtered$Genes)>=trimlen],
                    width(filtered$Genes)-trimlen,
                    width(filtered$Genes)),
               filtered$Genes[width(filtered$Genes)<trimlen])
    ind_genes <- c(which(width(filtered$Genes)>=trimlen), which(width(filtered$Genes)<trimlen))
    highly_expressed <- filtered$HE[ind_genes]
    depth_of_coverage <- filtered$Depth[ind_genes]
  } else if (trimside == "random"){
    genes_ok <- filtered$Genes[width(filtered$Genes)>=trimlen]
    start_lims <- width(genes_ok)-150+1
    start_inds <- lapply(start_lims,randomTriplet) %>% unlist()
    genes <- c(subseq(genes_ok,
                      start_inds,
                      start_inds+trimlen-1),
               filtered$Genes[width(filtered$Genes)<trimlen])
    ind_genes <- c(which(width(filtered$Genes)>=trimlen), which(width(filtered$Genes)<trimlen))
    highly_expressed <- filtered$HE[ind_genes]
    depth_of_coverage <- filtered$Depth[ind_genes]
  } else {
    stop("trimside must be set to \"start\",\"end\", or \"random\"")
  }

  # codon table
  codon_table <- codonTable(genes)

  # codon pair counts
  codon_pair_table <- getPairCounts(genes, genetic_code = genetic_code)

  if(!is.null(depth_of_coverage) & !all_metrics){
    gc <- sum(alphabetFrequency(genes)[,2:3]*depth_of_coverage)/sum(alphabetFrequency(genes)*depth_of_coverage)
    return(data.frame(CUBHE = getWeightedCUBHE(codon_table,
                                     highly_expressed,
                                     depth_of_coverage,
                                     genetic_code = genetic_code),
                      GC = gc,
                      GCdiv = abs(0.5-gc),
                      ConsistencyHE = getWeightedConsistency(codon_table,
                                                     highly_expressed,
                                                     depth_of_coverage,
                                                     genetic_code = genetic_code),
                      CUB = getWeightedCUB(codon_table,
                                   highly_expressed,
                                   depth_of_coverage,
                                   genetic_code = genetic_code),
                      CPB = NA,
                      FilteredSequences = filtered$Filtered,
                      nHE = sum(highly_expressed),
                      stringsAsFactors = FALSE))
  } else if (is.null(depth_of_coverage) & !all_metrics){
    gc <- sum(alphabetFrequency(genes)[,2:3])/sum(alphabetFrequency(genes))
    return(data.frame(CUBHE = getCUB(codon_table,
                                     highly_expressed,
                                     method = "MILC",
                                     genetic_code = genetic_code),
                      GC = gc,
                      GCdiv = abs(0.5-gc),
                      ConsistencyHE = getCUB(codon_table,
                                             highly_expressed,
                                             method = "consistency",
                                             genetic_code = genetic_code),
                      CUB = getCUB(codon_table,
                                   highly_expressed,
                                   method = "MILCgenomic",
                                   genetic_code = genetic_code),
                      CPB = getCPB(codon_pair_table),
                      FilteredSequences = filtered$Filtered,
                      nHE = sum(highly_expressed),
                      stringsAsFactors = FALSE))
  } else if(!is.null(depth_of_coverage) & all_metrics){
    gc <- sum(alphabetFrequency(genes)[,2:3]*depth_of_coverage)/sum(alphabetFrequency(genes)*depth_of_coverage)
    return(data.frame(MILC = getWeightedCUBHE(codon_table,
                                               highly_expressed,
                                               depth_of_coverage,
                                               genetic_code = genetic_code),
                      ENCprime = getWeightedCUBHE(codon_table,
                                              highly_expressed,
                                              depth_of_coverage,
                                              genetic_code = genetic_code,
                                              method="ENCprime"),
                      B = getWeightedCUBHE(codon_table,
                                                  highly_expressed,
                                                  depth_of_coverage,
                                                  genetic_code = genetic_code,
                                                  method="B"),
                      SCUO = getWeightedCUBHE(codon_table,
                                                  highly_expressed,
                                                  depth_of_coverage,
                                                  genetic_code = genetic_code,
                                                  method="SCUO"),
                      MCB = getWeightedCUBHE(codon_table,
                                                  highly_expressed,
                                                  depth_of_coverage,
                                                  genetic_code = genetic_code,
                                                  method="MCB"),
                      nHE = sum(highly_expressed),
                      stringsAsFactors = FALSE))
  } else if (is.null(depth_of_coverage) & all_metrics){
    gc <- sum(alphabetFrequency(genes)[,2:3])/sum(alphabetFrequency(genes))
    return(data.frame(MILC = getCUB(codon_table,
                                     highly_expressed,
                                     method = "MILC",
                                     genetic_code = genetic_code),
                      ENCprime = getCUB(codon_table,
                                                  highly_expressed,
                                                  genetic_code = genetic_code,
                                                  method="ENCprime"),
                      B = getCUB(codon_table,
                                           highly_expressed,
                                           genetic_code = genetic_code,
                                           method="B"),
                      SCUO = getCUB(codon_table,
                                              highly_expressed,
                                              genetic_code = genetic_code,
                                              method="SCUO"),
                      MCB = getCUB(codon_table,
                                             highly_expressed,
                                             genetic_code = genetic_code,
                                             method="MCB"),
                      nHE = sum(highly_expressed),
                      stringsAsFactors = FALSE))
  }
}



getCodonStatistics_i <- function(genes,
                                 highly_expressed,
                                 fragments = FALSE,
                                 depth_of_coverage = NULL,
                                 genetic_code = "11",
                                 n_le=100,
                                 trimlen = NA,
                                 trimside = "start",
                                 all_metrics=F){

  if(sum(highly_expressed) == 0){
    stop("No highly expressed genes?")
  }

  #Remove short sequences and sequences that are not multiples of 3
  if(!is.na(fragments)){
    if(fragments == TRUE){
      filtered <- filterSeq(genes = genes,
                            highly_expressed = highly_expressed,
                            length_threshold = 120,
                            depth_of_coverage = depth_of_coverage)
    } else {
      filtered <- filterSeq(genes = genes,
                            highly_expressed = highly_expressed,
                            depth_of_coverage = depth_of_coverage)
    }
  } else {
    filtered <- filterSeq(genes = genes,
                          highly_expressed = highly_expressed,
                          length_threshold = 1,
                          depth_of_coverage = depth_of_coverage)
  }

  if(is.na(trimlen)){
    genes <- filtered$Genes
    highly_expressed <- filtered$HE
    depth_of_coverage <- filtered$Depth
  } else if (trimside == "start"){
    genes <- c(subseq(filtered$Genes[width(filtered$Genes)>=trimlen], 1, trimlen),
               filtered$Genes[width(filtered$Genes)<trimlen])
    ind_genes <- c(which(width(filtered$Genes)>=trimlen), which(width(filtered$Genes)<trimlen))
    highly_expressed <- filtered$HE[ind_genes]
    depth_of_coverage <- filtered$Depth[ind_genes]
  } else if (trimside == "end"){
    genes <- c(subseq(filtered$Genes[width(filtered$Genes)>=trimlen],
                      width(filtered$Genes)-trimlen,
                      width(filtered$Genes)),
               filtered$Genes[width(filtered$Genes)<trimlen])
    ind_genes <- c(which(width(filtered$Genes)>=trimlen), which(width(filtered$Genes)<trimlen))
    highly_expressed <- filtered$HE[ind_genes]
    depth_of_coverage <- filtered$Depth[ind_genes]
  } else if (trimside == "random"){
      genes_ok <- filtered$Genes[width(filtered$Genes)>=trimlen]
      start_lims <- width(genes_ok)-trimlen+1
      start_inds <- lapply(start_lims,randomTriplet) %>% unlist()
      genes <- c(subseq(genes_ok,
                        start_inds,
                        start_inds+trimlen-1),
                 filtered$Genes[width(filtered$Genes)<trimlen])
      ind_genes <- c(which(width(filtered$Genes)>=trimlen), which(width(filtered$Genes)<trimlen))
      highly_expressed <- filtered$HE[ind_genes]
      depth_of_coverage <- filtered$Depth[ind_genes]
  } else {
    stop("trimside must be set to \"start\",\"end\", or \"random\"")
  }

  # codon table
  codon_table <- codonTable(genes)

  # codon pair counts
  codon_pair_table <- getPairCounts(genes, genetic_code = genetic_code)

  if(!is.null(depth_of_coverage) & !all_metrics){
    gc <- sum(alphabetFrequency(genes)[,2:3]*depth_of_coverage)/sum(alphabetFrequency(genes)*depth_of_coverage)
    cubi <- getWeightedCUBi(genes,highly_expressed,depth_of_coverage,n_le=n_le)
    return(data.frame(CUBHE = cubi[2],
                      GC = gc,
                      GCdiv = abs(0.5-gc),
                      ConsistencyHE = getWeightedConsistency(codon_table,
                                                             highly_expressed,
                                                             depth_of_coverage,
                                                             genetic_code = genetic_code),
                      CUB = cubi[1],
                      CPB = NA,
                      FilteredSequences = filtered$Filtered,
                      nHE = sum(highly_expressed),
                      stringsAsFactors = FALSE))
  } else if(is.null(depth_of_coverage) & !all_metrics){
    cubi <- CUBi(genes,highly_expressed,n_le=n_le)
    gc <- sum(alphabetFrequency(genes)[,2:3])/sum(alphabetFrequency(genes))
    return(data.frame(CUBHE = cubi[2],
                      GC = gc,
                      GCdiv = abs(0.5-gc),
                      ConsistencyHE = getCUB(codon_table,
                                             highly_expressed,
                                             method = "consistency",
                                             genetic_code = genetic_code),
                      CUB = cubi[1],
                      CPB = NA,
                      FilteredSequences = filtered$Filtered,
                      nHE = sum(highly_expressed),
                      stringsAsFactors = FALSE))
  } else if(!is.null(depth_of_coverage) & all_metrics){
    gc <- sum(alphabetFrequency(genes)[,2:3]*depth_of_coverage)/sum(alphabetFrequency(genes)*depth_of_coverage)
    cubi <- getWeightedCUBi(genes,highly_expressed,depth_of_coverage,n_le=n_le)
    return(data.frame(MILC = cubi[2],
                      ENCprime = getWeightedCUBi(genes,highly_expressed,depth_of_coverage,n_le=n_le,
                                                 method="ENCprime")[2],
                      B = getWeightedCUBi(genes,highly_expressed,depth_of_coverage,n_le=n_le,
                                                 method="B")[2],
                      SCUO = getWeightedCUBi(genes,highly_expressed,depth_of_coverage,n_le=n_le,
                                                 method="SCUO")[2],
                      MCB = getWeightedCUBi(genes,highly_expressed,depth_of_coverage,n_le=n_le,
                                                 method="MCB")[2],
                      nHE = sum(highly_expressed),
                      stringsAsFactors = FALSE))
  } else if(is.null(depth_of_coverage) & all_metrics){
    cubi <- CUBi(genes,highly_expressed,n_le=n_le)
    gc <- sum(alphabetFrequency(genes)[,2:3])/sum(alphabetFrequency(genes))
    return(data.frame(MILC = cubi[2],
                      ENCprime = CUBi(genes,highly_expressed,n_le=n_le,
                                      method="ENCprime")[2],
                      B = CUBi(genes,highly_expressed,n_le=n_le,
                                      method="B")[2],
                      SCUO = CUBi(genes,highly_expressed,n_le=n_le,
                                      method="SCUO")[2],
                      MCB = CUBi(genes,highly_expressed,n_le=n_le,
                                      method="MCB")[2],
                      nHE = sum(highly_expressed),
                      stringsAsFactors = FALSE))
  }
}
