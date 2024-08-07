

source('./startup.R')


CDROM_func <- read.csv("./Dup_Functionalizations.tsv", sep="") %>%
  select(dup_1 = Dup1, dup_2 = Dup2, anc = Ancestor, CDROM_func = func) %>%
  na.omit()


########
species_pairs <- c('AN_YA', "ME_YA", 'MO_VI', 'PE_PS', 'PS_PE', "VI_MO", 'WI_PS', "YA_ME")



get_dups_for_species_pair <- function(species_1, species_2) {
  orthogroups <- read.delim2("./CLOUD/Orthogroups.tsv", na.strings = '')
  orig_exp <- read.csv("./CLOUD/Expression_Data.tsv", sep="")
  
  
  # keep only orthogroups that have trees so that I can get branch lengths for cloud input 
  og_with_tree <- orthogroups %>%
    mutate(n_na = rowSums(is.na(across(everything())))) %>%
    filter(n_na <= 5)
  
  dups_anc <- read.csv("./CLOUD/Dup_Pairs_Ancestral.tsv", sep="")
  dups_anc <- dups_anc %>%
    filter(Orthogroup %in% og_with_tree$Orthogroup)
  
  
  # format eight species into duplicate + ancestral species pairs 
  dups_anc <- dups_anc %>%
    mutate(dup_species = gsub("\\d|YOgn", "", dup_1),
           anc_species = gsub("\\d|YOgn", "", ancestral_copy)) %>%
    filter((str_detect(dup_species, species_1) & str_detect(anc_species, species_2)) |
             (str_detect(dup_species, species_2) & str_detect(anc_species, species_1)))
  
  # extract the branch lengths (tpc and tpca for cloud input)
  dups_anc$TPC <- NA
  dups_anc$TPCA <- NA
  
  p = 0
  for(i in 1:nrow(dups_anc)) {
    row <- dups_anc[i,]
    
    OG <- row$Orthogroup
    dup1 <- row$dup_1
    dup2 <- row$dup_2
    anc <- row$ancestral_copy
    
    
    tree <- read.tree(paste0('./CLOUD/Gene_Trees/', OG, '_tree.txt'))
    
    dup1_length <- setNames(tree$edge.length[sapply(1:length(tree$tip.label),function(x,y) which(y==x),y=tree$edge[,2])],tree$tip.label)[dup1]
    dup2_length <- setNames(tree$edge.length[sapply(1:length(tree$tip.label),function(x,y) which(y==x),y=tree$edge[,2])],tree$tip.label)[dup2]
    longer_dup <- ifelse(dup1_length > dup2_length, names(dup1_length), names(dup2_length))
    
    #if (cophenetic(tree)[dup1, dup2] != dup1_length + dup2_length){print('false')}
    
    tpc <- cophenetic(tree)[dup1, dup2]
    tpca <- cophenetic(tree)[longer_dup, anc]
    
    dups_anc$TPC[i] <- tpc
    dups_anc$TPCA[i] <- tpca
    
    if(tpc > tpca) {p = p+1}
    
  }
  
  # have to remove completely identical duplicates bc tpc and maybe even tpca will be 0 
  dups_anc <- dups_anc %>%
    filter(TPC > 0 & TPCA > 0)
  
  # run cloud on each species pair 
  func_output <- data.frame()
  
  # get duplicate+anc species pairs for the given species pair
  dup_input <- dups_anc %>%
    #filter(species_pair == pair) %>%
    select(dup_1,TPC,TPCA,dup_2,ancestral_copy)
  
  exp <- orig_exp
  orig_dup_input <- dup_input
  
  # merge dups with expression data
  colnames(exp) <- paste0(colnames(orig_exp),'_d1')
  dup_input <- merge(dup_input, exp, by.x='dup_1', by.y='YOgnID_d1')
  
  colnames(exp) <- paste0(colnames(orig_exp),'_d2')
  dup_input <- merge(dup_input, exp, by.x='dup_2', by.y='YOgnID_d2')
  
  
  # format duplicate gene cloud input
  colnames(exp) <- paste0(colnames(orig_exp),'_anc')
  dup_input <- merge(dup_input, exp, by.x='ancestral_copy', by.y='YOgnID_anc')
  
  dup_input <- dup_input %>%
    select(dup_1, dup_2, ancestral_copy)
  
  return(dup_input)
}

results <- data.frame()

for (pair in species_pairs) {
  
  files <- list.files(path = "./CLOUD/output_CLOUD/", pattern = pair, full.names = TRUE)
  
  species_pair <- strsplit(pair, "_")[[1]]
  species_1 <- species_pair[1]
  species_2 <- species_pair[2]
  
  CLOUD_func <- lapply(files, read.table, header = FALSE, stringsAsFactors = FALSE) %>% bind_cols()
  
  if (ncol(CLOUD_func) == 3) {CLOUD_func$class_4 <- NA}
  
  colnames(CLOUD_func) <- c('class_1', 'class_2', 'class_3', 'class_4')
  CLOUD_func <- CLOUD_func[-1,]
  
  dups <- get_dups_for_species_pair(species_1, species_2)
  
  CLOUD_func <- cbind(dups, CLOUD_func)
  
  results <- rbind(results, CLOUD_func)
  
  
}

########

CLOUD_func <- results %>%
  distinct(dup_1, dup_2, ancestral_copy, .keep_all = T)

CLOUD_CDROM_func <- merge(CLOUD_func, CDROM_func, by = c('dup_1', 'dup_2'))


# keep only the duplicate pairs with CLOUD funcs the same for all 3 runs 
max_identical <- function(...) {max(table(unlist(list(...))))}
CLOUD_CDROM_func <- CLOUD_CDROM_func %>% 
  rowwise() %>% 
  select(-class_4) %>%
  mutate(max_identical = max_identical(c_across(class_1:class_3))) %>%
  mutate_all(~str_replace_all(., "parent", "_dup1")) %>%
  mutate_all(~str_replace_all(., "child", "_dup2"))

table(CLOUD_CDROM_func$max_identical)

all_CLOUD_same <- CLOUD_CDROM_func %>%
  filter(max_identical == 3)

# compare CLOUD and CDROM funcs
table(all_CLOUD_same$class_1, all_CLOUD_same$CDROM_func)

# number of duplicate pairs classified into the same funcs by CDROM and CLOUD 
sum(all_CLOUD_same$class_1 == all_CLOUD_same$CDROM_func)

# percentage:
(sum(all_CLOUD_same$class_1 == all_CLOUD_same$CDROM_func) / nrow(all_CLOUD_same)) * 100

# format and write to file CLOUD+CDROM funcs dataframe
all_CLOUD_same <- all_CLOUD_same %>%
  mutate(CLOUD_func = class_1) %>%
  select(dup_1, dup_2, anc, CDROM_func, CLOUD_func)

write.table(all_CLOUD_same, file = 'CDROM_CLOUD_funcs.tsv')


# add pseudogenization (classified in script 2_OrthoFinder_Output_Cleanup)
pseudo <- read.csv("./Dup_Functionalizations.tsv", sep="") %>% 
  select(dup_1 = Dup1, dup_2 = Dup2, anc = Ancestor, pseudo = pseudo)

CLOUD_CDROM_pseudo <- left_join(pseudo, all_CLOUD_same, by = c('dup_1', 'dup_2', 'anc'))

CLOUD_CDROM_pseudo <- CLOUD_CDROM_pseudo %>%
  mutate(CDROM_pseudo = coalesce(CDROM_func, pseudo),
         CLOUD_pseudo = coalesce(CLOUD_func, pseudo)) %>%
  filter(!is.na(CDROM_pseudo))

write.table(CLOUD_CDROM_pseudo, file = 'CDROM_CLOUD_pseudo_funcs.tsv')


