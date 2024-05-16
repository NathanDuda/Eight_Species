
source('./startup.R')

dups <- read.csv("./Duplicate_Pairs.tsv", sep="")
longest_transcript <- read.csv("./longest_transcript.tsv", sep="") %>%
  select(YOgn, longest_ORF)



colnames(longest_transcript) <- c('dup_1', 'dup_1_seq')
dup_seqs <- merge(dups, longest_transcript, by = 'dup_1')

colnames(longest_transcript) <- c('dup_2', 'dup_2_seq')
dup_seqs <- merge(dup_seqs, longest_transcript, by = 'dup_2')


# 

calculate_mutations <- function(seq1, seq2) {
  # Function to calculate mutations
  nonsynonymous <- 0
  synonymous <- 0
  nonsense1 <- 0
  nonsense2 <- 0
  
  nonsynonymous_positions <- c()
  synonymous_positions <- c()
  nonsense1_positions <- c()
  nonsense2_positions <- c()
  
  
  if (nchar(seq1) %% 3 != 0) {stop('length of seq1 is not a a multiple of 3')}
  if (nchar(seq2) %% 3 != 0) {stop('length of seq2 is not a a multiple of 3')}
  
  # Check length of sequences
  if (nchar(seq1) != nchar(seq2)) {
    stop("Sequences are not aligned properly")
  }
  
  for (i in 1:nchar(seq1)) {
    position = (i * 3) - 3
    codon1 <- substr(seq1, position + 1, position + 3)
    codon2 <- substr(seq2, position + 1, position + 3)
    
    # Check if codons are complete
    if (nchar(codon1) == 3 && nchar(codon2) == 3) {
      # Check for synonymous mutation
      if (codon1 != codon2) {
        if (translate_codon(codon1) == translate_codon(codon2)) {
          synonymous <- synonymous + 1
          synonymous_positions <- c(synonymous_positions, position)
        }
        
        stop_codons <- c("TAA", "TAG", "TGA")
        if (codon1 %in%  stop_codons) {
          nonsense1 <- nonsense1 + 1
          nonsense1_positions <- c(nonsense1_positions, position)
        }
        if (codon2 %in%  stop_codons) {
          nonsense2 <- nonsense2 + 1
          nonsense2_positions <- c(nonsense2_positions, position)
        }
        
        if ((translate_codon(codon1) != translate_codon(codon2))) { # nonsense mutations are counted in nonsyn too
          nonsynonymous <- nonsynonymous + 1
          nonsynonymous_positions <- c(nonsynonymous_positions, position)
        }
        #if ((translate_codon(codon1) != translate_codon(codon2)) & (!codon2 %in% stop_codons)) {
        #  nonsynonymous2 <- nonsynonymous2 + 1
        #  nonsynonymous2_positions <- c(nonsynonymous2_positions, position)
        #}
        
        
      }
    }
  }
  
  
  return(list(nonsynonymous = nonsynonymous, 
              synonymous = synonymous, 
              nonsense1 = nonsense1,
              nonsense2 = nonsense2,
              
              nonsynonymous_positions = nonsynonymous_positions, 
              synonymous_positions = synonymous_positions, 
              nonsense1_positions = nonsense1_positions,
              nonsense2_positions = nonsense2_positions))
}


translate_codon <- function(codon) {
  # Function to translate codon to amino acid
  if (grepl("-", codon)) {return("-")}
  
  codon_table <- list(
    "TTT" = "F", "TTC" = "F", "TTA" = "L", "TTG" = "L",
    "TCT" = "S", "TCC" = "S", "TCA" = "S", "TCG" = "S",
    "TAT" = "Y", "TAC" = "Y", "TAA" = "*", "TAG" = "*",
    "TGT" = "C", "TGC" = "C", "TGA" = "*", "TGG" = "W",
    "CTT" = "L", "CTC" = "L", "CTA" = "L", "CTG" = "L",
    "CCT" = "P", "CCC" = "P", "CCA" = "P", "CCG" = "P",
    "CAT" = "H", "CAC" = "H", "CAA" = "Q", "CAG" = "Q",
    "CGT" = "R", "CGC" = "R", "CGA" = "R", "CGG" = "R",
    "ATT" = "I", "ATC" = "I", "ATA" = "I", "ATG" = "M",
    "ACT" = "T", "ACC" = "T", "ACA" = "T", "ACG" = "T",
    "AAT" = "N", "AAC" = "N", "AAA" = "K", "AAG" = "K",
    "AGT" = "S", "AGC" = "S", "AGA" = "R", "AGG" = "R",
    "GTT" = "V", "GTC" = "V", "GTA" = "V", "GTG" = "V",
    "GCT" = "A", "GCC" = "A", "GCA" = "A", "GCG" = "A",
    "GAT" = "D", "GAC" = "D", "GAA" = "E", "GAG" = "E",
    "GGT" = "G", "GGC" = "G", "GGA" = "G", "GGG" = "G"  
  )
  
  return(codon_table[[codon]])
}



#


results <- data.frame() 
for (row in 1:nrow(dup_seqs)) {
  
  pair <- dup_seqs[row,]
  
  seq1 <- pair$dup_1_seq
  seq2 <- pair$dup_2_seq
  seq1_name <- pair$dup_1
  seq2_name <- pair$dup_2
  
  # align seqs
  alignment <- pwalign::pairwiseAlignment(pattern = seq1, subject = seq2, substitutionMatrix = "BLOSUM62")
  
  seq1 <- alignment@pattern
  seq2 <- alignment@subject
  
  
  # make sure aligned sequences are length multiple of 3 (add '-' to end when not)
  num_dashes <- (3 - (nchar(seq1) %% 3)) %% 3
  if (num_dashes == 1) {seq1 <- paste0(as.character(seq1), '-')}
  if (num_dashes == 2) {seq1 <- paste0(as.character(seq1), '--')}
  
  num_dashes <- (3 - (nchar(seq2) %% 3)) %% 3
  if (num_dashes == 1) {seq2 <- paste0(as.character(seq2), '-')}
  if (num_dashes == 2) {seq2 <- paste0(as.character(seq2), '--')}
  
  
  # find mutations
  mutations <- calculate_mutations(seq1, seq2)
  
  
  # format results
  df <- data.frame(
    type = c("Nonsynonymous", "Synonymous", "Nonsense1", "Nonsense2"),
    count = c(mutations$nonsynonymous, mutations$synonymous, mutations$nonsense1, mutations$nonsense2),
    positions = I(list(mutations$nonsynonymous_positions, mutations$synonymous_positions, mutations$nonsense1_positions, mutations$nonsense2_positions)),
    seq1 = seq1_name,
    seq2 = seq2_name,
    aligned_seq_length = nchar(seq1)
  )
  
  # combine results
  results <- rbind(results, df)
  
  print(paste0(row / nrow((dup_seqs)) * 100, ' %'))
  
}

results$positions <- sapply(results$positions, paste, collapse = ",")
results <- as.data.frame(results)
write.table(results, file = './duppairs_nonsense_dnds_mutations.tsv')





nonsense <- read.csv2("./duppairs_nonsense_dnds_mutations.tsv", sep="")


nonsense <- nonsense %>%
  group_by(seq1) %>%
  mutate(has_nonsense = any((type == 'Nonsense1' | type == 'Nonsense2') & count != 0)) %>%
  ungroup() %>%
  select(type, count, seq1, seq2, has_nonsense)

nonsense %>%
  #pivot_wider(names_from = 'type', values_from = 'count')
  ggplot(aes(x = type, y = count, color = has_nonsense)) +
    geom_boxplot()


t <- nonsense %>%
  filter(type == 'Nonsense1' | type == 'Nonsense2') %>%
  filter(count != 0)












