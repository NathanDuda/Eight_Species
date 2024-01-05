


source('./startup.R')

# read in duplicate pairs with their ancestral copy 
dups <- read.csv("./Dup_Pairs_Ancestral.tsv", sep="")

# read in nucleotide and protein sequences for all genes
gn_seqs <- read.csv("./longest_transcript.tsv", sep="")
gn_seqs <- gn_seqs[c('YOgn','longest_ORF','prot')]

# read in all one to two orthogroups with their genes per species 
orthogroups <- read.csv("./Dup_Pair_Orthologs.tsv", sep="")

# make a list of all species names 
species_list <- c('dana','dmel','dmoj','dper','dpse','dvir','dwil','dyak')



### one to two orthogroups:

# write nucleotide fasta files 
for (row_num in nrow(orthogroups):1) {
  orthogroup <- orthogroups[row_num,]
  orthogroup_name <- orthogroup$Orthogroup
  orthogroup <- orthogroup %>% 
    separate((which(species_list == orthogroup$duplicate_pair_species) + 1), sep = ", ", into = c('a','b')) %>%
    select(-Orthogroup, -duplicate_pair_species)
  for (gene in orthogroup) {
    if (nchar(gene) != 0) {
      gene_seq <- gn_seqs[which(gn_seqs$YOgn == gene),2]
      output_file <- paste0('./Grouped_Fastas/One_to_Two_Orthogroups/Nucleotide_Fastas/',orthogroup_name,'.fa')  
      cat(">", gene, '\n', gene_seq, "\n", file = output_file, append = T, sep = '')
    }
  }
}

# write protein fasta files 
for (row_num in nrow(orthogroups):1) {
  orthogroup <- orthogroups[row_num,]
  orthogroup_name <- orthogroup$Orthogroup
  orthogroup <- orthogroup %>% 
    separate((which(species_list == orthogroup$duplicate_pair_species) + 1), sep = ", ", into = c('a','b')) %>%
    select(-Orthogroup, -duplicate_pair_species)
  for (gene in orthogroup) {
    if (nchar(gene) != 0) {
      gene_seq <- gn_seqs[which(gn_seqs$YOgn == gene),3]
      output_file <- paste0('./Grouped_Fastas/One_to_Two_Orthogroups/Protein_Fastas/',orthogroup_name,'.fa')  
      cat(">", gene, '\n', gene_seq, "\n", file = output_file, append = T, sep = '')
    }
  }
}



### one to one orthogroups:
one_to_ones <- read.csv("./One_to_Ones.tsv", sep="")

# write nucleotide fasta files
for (row_num in nrow(one_to_ones):1) {
  orthogroup <- one_to_ones[row_num,]
  orthogroup_name <- orthogroup$Orthogroup
  orthogroup <- orthogroup %>% select(-Orthogroup)
  for (gene in orthogroup) {
    if (nchar(gene) != 0) {
      gene_seq <- gn_seqs[which(gn_seqs$YOgn == gene),2]
      output_file <- paste0('./Grouped_Fastas/One_to_One_Orthogroups/Nucleotide_Fastas/',orthogroup_name,'.fa')  
      cat(">", gene, '\n', gene_seq, "\n", file = output_file, append = T, sep = '')
    }
  }
}

# write protein fasta files 
for (row_num in nrow(one_to_ones):1) {
  orthogroup <- one_to_ones[row_num,]
  orthogroup_name <- orthogroup$Orthogroup
  orthogroup <- orthogroup %>% select(-Orthogroup)
  for (gene in orthogroup) {
    if (nchar(gene) != 0) {
      gene_seq <- gn_seqs[which(gn_seqs$YOgn == gene),3]
      output_file <- paste0('./Grouped_Fastas/One_to_One_Orthogroups/Protein_Fastas/',orthogroup_name,'.fa')  
      cat(">", gene, '\n', gene_seq, "\n", file = output_file, append = T, sep = '')
    }
  }
}



### duplicate and ancestral pairings:

# write nucleotide and protein fasta files for dup_1 with anc, dup_2 with anc, and dup_1 with dup_2 pairings
for (row_num in nrow(dups):1) {
  # get current set of dup_1, dup_2, and anc 
  dup1_dup2_anc <- one_to_ones[row_num,]
  
  # combine the pair gene names to get the output fasta file name 
  dup1_anc_name <- paste0(dup1_dup2_anc$dup_1,'_',dup1_dup2_anc$ancestral_copy)
  dup2_anc_name <- paste0(dup1_dup2_anc$dup_2,'_',dup1_dup2_anc$ancestral_copy)
  dup1_dup2_name <- paste0(dup1_dup2_anc$dup_1,'_',dup1_dup2_anc$dup_2)
  
  dup1_dup2_anc <- dup1_dup2_anc %>% select(-Orthogroup)
  
  # extract the gene name for each gene copy 
  dup_1 <- dup1_dup2_anc$dup_1
  dup_2 <- dup1_dup2_anc$dup_2
  anc <- dup1_dup2_anc$ancestral_copy
  
  # extract the nucleotide sequences for each gene copy 
  dup1_nuc <- gn_seqs[which(gn_seqs$YOgn == dup_1),2]
  dup2_nuc <- gn_seqs[which(gn_seqs$YOgn == dup_2),2]
  anc_nuc <- gn_seqs[which(gn_seqs$YOgn == anc),2]
  
  # extract the protein sequences for each gene copy 
  dup1_prot <- gn_seqs[which(gn_seqs$YOgn == dup_1),3]
  dup2_prot <- gn_seqs[which(gn_seqs$YOgn == dup_2),3]
  anc_prot <- gn_seqs[which(gn_seqs$YOgn == anc),3]
  
  
  # write the nucleotide fastas for each pairing
  output_file <- paste0('./Grouped_Fastas/Dup1_Anc/Nucleotide_Fastas/',dup1_anc_name,'.fa')  
  cat(">", dup_1, '\n', dup1_nuc, "\n", file = output_file, append = T, sep = '')
  cat(">", anc, '\n', anc_nuc, "\n", file = output_file, append = T, sep = '')
  
  output_file <- paste0('./Grouped_Fastas/Dup2_Anc/Nucleotide_Fastas/',dup2_anc_name,'.fa')  
  cat(">", dup_2, '\n', dup2_nuc, "\n", file = output_file, append = T, sep = '')
  cat(">", anc, '\n', anc_nuc, "\n", file = output_file, append = T, sep = '')
  
  output_file <- paste0('./Grouped_Fastas/Dup1_Dup2/Nucleotide_Fastas/',dup1_dup2_name,'.fa')  
  cat(">", dup_1, '\n', dup1_nuc, "\n", file = output_file, append = T, sep = '')
  cat(">", anc, '\n', dup2_nuc, "\n", file = output_file, append = T, sep = '')
  
  
  # write the protein fastas for each pairing
  output_file <- paste0('./Grouped_Fastas/Dup1_Anc/Protein_Fastas/',dup1_anc_name,'.fa')  
  cat(">", dup_1, '\n', dup1_prot, "\n", file = output_file, append = T, sep = '')
  cat(">", anc, '\n', anc_prot, "\n", file = output_file, append = T, sep = '')
  
  output_file <- paste0('./Grouped_Fastas/Dup2_Anc/Protein_Fastas/',dup2_anc_name,'.fa')  
  cat(">", dup_2, '\n', dup2_prot, "\n", file = output_file, append = T, sep = '')
  cat(">", anc, '\n', anc_prot, "\n", file = output_file, append = T, sep = '')
  
  output_file <- paste0('./Grouped_Fastas/Dup1_Dup2/Protein_Fastas/',dup1_dup2_name,'.fa')  
  cat(">", dup_1, '\n', dup1_prot, "\n", file = output_file, append = T, sep = '')
  cat(">", anc, '\n', dup2_prot, "\n", file = output_file, append = T, sep = '')
}


