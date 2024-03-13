
# this script also cleans up the raw expression data 
# and calculates tau for every expressed gene


source('./startup.R')

# OrthoFinder was run on with default parameters the protein fastas:
# OrthoFinder/orthofinder -f ./Eight_Species/Protein_Fastas/ -o ./Eight_Species/OrthoFinder_Output/


longest_transcript <- read.csv("./longest_transcript.tsv", sep="")



### Orthologs:
orthogroups <- read.delim2("./OrthoFinder_Output/Results_Jan01/Orthogroups/Orthogroups.tsv")

# get the one to ones 
one_to_ones <- orthogroups %>%
  filter_all(all_vars(!grepl(",", .))) %>% 
  mutate_all(~ifelse(. == "", NA, .)) %>% # change empty cells to NA 
  na.omit()

# write the one to ones to file
write.table(one_to_ones,file='One_to_Ones.tsv')


# import the number of genes each species has in each orthogroup 
orthogroup_gene_count <- read.delim("./OrthoFinder_Output/Results_Jan01/Orthogroups/Orthogroups.GeneCount.tsv")

# get the two to one to ones to zeros 
two_to_ones <- orthogroup_gene_count %>%
  filter(if_all(2:9, ~. <= 2)) %>%              # keep only pairs 
  filter(rowSums(select(., 2:9) == 2) == 1) %>% # make sure there is only one species with 2 copies 
  filter(Total > 2)                             # remove duplicates without any orthologs 
  
  
# add a column with the name of the species with the duplication
two_to_ones$duplicate_pair_species <- 
  apply(two_to_ones[, 2:9], 1, function(x) {
  col_index <- which(x == 2)
  return(colnames(two_to_ones[, 2:9])[col_index])})

two_to_ones <- two_to_ones[c('Orthogroup','duplicate_pair_species')]
two_to_ones$duplicate_pair_species <- gsub('_prot','',two_to_ones$duplicate_pair_species)

# merge back with the gene names 
two_to_ones <- merge(orthogroups,two_to_ones,by='Orthogroup')

write.table(two_to_ones, file='Dup_Pair_Orthologs.tsv')



### Paralogs: 

dups <- two_to_ones %>%
  rowwise() %>%
  
  # extract the duplicate pair genes
  mutate(duplicate_pair = toString(c_across(2:9)[grep(",", c_across(2:9))])) %>%
  select(Orthogroup, duplicate_pair) %>%
  separate(duplicate_pair, into = c("dup_1", "dup_2"), sep = ", ") %>%
  
  # ordering so that dup_1 is always > dup_2 for easy error catching  
  mutate(temp_column = dup_1) %>%
  mutate(dup_1 = case_when(dup_1 < dup_2 ~ dup_1, dup_1 > dup_2 ~ dup_2)) %>%
  mutate(dup_2 = case_when(dup_1 < dup_2 ~ dup_2, dup_1 > dup_2 ~ temp_column)) %>%
  select(-temp_column)

# keep only duplicate pairs with expression 
species_names <- c('dana','dmel','dmoj','dper','dpse','dvir','dwil','dyak')

all_expression <- data.frame()
for (species in species_names){
  
  expression <- read.delim(paste0('./Raw_Data/YO_Expression/GSE99574_HiSAT2_',species,'.nrc.YO.txt'))
  expression <- expression %>% select(-jaccard, -FBgnID)
  
  # remove the unique species names in dmel 
  if (species == 'dmel') {
    expression <- expression %>%
      rename_all(~ ifelse(startsWith(., "w1118"), paste0(., "1"), .)) %>% # prevent duplicate column names
      rename_all(~ gsub(paste0('w1118_'), "", .)) %>%
      rename_all(~ gsub(paste0('oreR_'), "", .))
    }
  
  expression <- expression %>% 
    rename_all(~ gsub(paste0(species,"_"), "", .)) %>%
    mutate(across(2:ncol(expression), as.numeric)) 
  
  # get average expression value for replicates
  tissue_names <- c('f_ac','f_dg','f_go','f_hd','f_re','f_tx','f_wb',
                    'm_ac','m_dg','m_go','m_hd','m_re','m_tx','m_wb')
  
  for (tissue in tissue_names) {
    expression <- expression %>%
      mutate("{tissue}" := rowMeans(select(., starts_with(tissue)), na.rm = TRUE)) # combine replicates by avg
      #mutate("{tissue}" := rowMeans(select(., starts_with(tissue) & ends_with('R1')), na.rm = TRUE)) # use only one replicate for test
  }
  
  # filter out genes with no expression in all tissues 
  expressed_genes <- expression %>%
    select(YOgnID,all_of(tissue_names)) %>%
    filter(rowSums(select(., 2:15)) > 0)
  
  all_expression <- rbind(all_expression,expressed_genes)
}

# make expression values lower than 1 into 0 
all_expression[all_expression < 1] <- 0 
# re filter genes without any expression 
all_expression <- all_expression %>%
  filter(rowSums(select(., 2:15)) > 0)

# write expressed genes to file
write.table(all_expression,'Expression_Data.tsv')


# keep only expressed duplicates 
colnames(all_expression)[1] <- 'dup_1'
dups_expressed <- merge(dups,all_expression,by='dup_1')

colnames(all_expression)[1] <- 'dup_2'
dups_expressed <- merge(dups_expressed,all_expression,by='dup_2')

# write duplicate pairs to file
dups <- dups_expressed[c('Orthogroup','dup_1','dup_2')]
write.table(dups, file = 'Duplicate_Pairs.tsv')



# calculating tau for every expressed gene

# remotes::install_github("roonysgalbi/tispec")
library(tispec)

colnames(all_expression)[1] <- 'YOgnID'
rownames(all_expression) <- all_expression$YOgnID
all_expression <- all_expression %>% select(-YOgnID)

tau <- calcTau(all_expression)

tau$YOgnID <- rownames(tau)
tau <- tau[c('YOgnID','tau')]

# write tau to file
write.table(tau, file= 'Tau.tsv')



# sex bias 

library(DESeq2)

exp_deseq <- all_expression 

col_data <- data.frame(sex = as.factor(c(rep("female", 7), rep("male", 7))))
rownames(col_data) <- colnames(exp_deseq)


exp_deseq <- round(exp_deseq, digits = 0)

dds <- DESeqDataSetFromMatrix(countData = exp_deseq, colData = col_data, design = ~ sex)
dds <- DESeq(dds)
res <- results(dds)

male_biased <- rownames(res)[which(res$padj < 0.05 & res$log2FoldChange > 0)]
female_biased <- rownames(res)[which(res$padj < 0.05 & res$log2FoldChange < 0)]

exp_deseq$bias <- "neutral"
exp_deseq$bias[rownames(exp_deseq) %in% male_biased] <- "male_biased"
exp_deseq$bias[rownames(exp_deseq) %in% female_biased] <- "female_biased"

sex_bias <- exp_deseq %>%
  mutate(YOgnID = rownames(.)) %>%
  select(YOgnID, bias)

write.table(sex_bias, file = './sex_bias.tsv')
