
countCodonPairs <- function(codon_vec){
  codon_pairs <- data.frame(C1=codon_vec[1:(length(codon_vec) - 1)],
                            C2=codon_vec[2:(length(codon_vec))],
                            stringsAsFactors = F)
  codon_pairs$P <- paste(codon_pairs$C1, codon_pairs$C2, sep = "-")
  codon_pairs <- codon_pairs %>% group_by(P) %>% dplyr::count()
  codon_pairs$Codon1 <- codon_pairs$P %>% gsub(pattern = "-.*", replace = "")
  codon_pairs$Codon2 <- codon_pairs$P %>% gsub(pattern = ".*-", replace = "")
  return(codon_pairs)
}

getAA <- function(codon_pairs, genetic_code = "11"){
  codon_pairs$AA1 <- codon_pairs$Codon1 %>% DNAStringSet() %>%
    translate(genetic.code = getGeneticCode(id_or_name2=genetic_code)) %>%
    as.character()
  codon_pairs$AA2 <- codon_pairs$Codon2 %>% DNAStringSet() %>%
    translate(genetic.code = getGeneticCode(id_or_name2=genetic_code)) %>%
    as.character()
  codon_pairs$PAA <- paste(codon_pairs$AA1, codon_pairs$AA2, sep="-")
  return(codon_pairs)
}

getPairCounts <- function(genes, genetic_code = "11"){
  genome <- unlist(genes)

  #Group Nucleotides into codons
  codon_vec <- genome %>% codons() %>% as.character()
  #Get AA

  aa_vec <- translate(genome,
                      genetic.code = getGeneticCode(id_or_name2=genetic_code)) %>%
    as.character() %>% strsplit(split="") %>% unlist()

  #Table of codon pair counts w/ columns for Codon and AA ids
  count_tbl <- codon_vec %>% countCodonPairs() %>%
    getAA(.,genetic_code = genetic_code) %>% as.data.frame(stringsAsFactors = F)
  count_tbl$n <- as.numeric(count_tbl$n)

  #Single codon couns
  codon_tbl <- table(codon_vec)

  #Single aa counts
  aa_tbl <- table(aa_vec)

  #aa pair counts
  aa_pair_tbl <- count_tbl %>% subset(select = c(PAA,n)) %>%
    group_by(PAA) %>% summarise(PAA_n = sum(n)) %>%
    as.data.frame(stringsAsFactors = F)
  rownames(aa_pair_tbl) <- aa_pair_tbl$PAA

  #Attach other counts to codon pair counts
  count_tbl$Codon1_n <- codon_tbl[count_tbl$Codon1] %>% as.numeric()
  count_tbl$Codon2_n <- codon_tbl[count_tbl$Codon2] %>% as.numeric()
  count_tbl$AA1_n <- aa_tbl[count_tbl$AA1] %>% as.numeric()
  count_tbl$AA2_n <- aa_tbl[count_tbl$AA2] %>% as.numeric()
  count_tbl$PAA_n <- aa_pair_tbl[count_tbl$PAA,"PAA_n"] %>% as.numeric()
  return(count_tbl)
}

getCPB <- function(count_tbl){
  count_tbl <- count_tbl %>%
    mutate(CPS = log(n / (PAA_n * (Codon1_n * Codon2_n) / (AA1_n * AA2_n))))

  CPB <- sum(count_tbl$CPS) / (nrow(count_tbl) - 1)
  return(CPB)
}
