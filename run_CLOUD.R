#!/usr/bin/env Rscript


# on cluster:

#source('./startup.R')

pair <- commandArgs(trailingOnly = TRUE)




species_pair <- strsplit(pair, "_")[[1]]
species_1 <- species_pair[1]
species_2 <- species_pair[2]


source("./CLOUD/CLOUD.r")


library(ape)
library(dplyr)
library(tensorflow)
# install_tensorflow()
library(keras) 
library(mvtnorm)
library(stringr)


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
  
  if (cophenetic(tree)[dup1, dup2] != dup1_length + dup2_length){print('false')}
  
  tpc <- cophenetic(tree)[dup1, dup2]
  tpca <- cophenetic(tree)[longer_dup, anc]
  
  dups_anc$TPC[i] <- tpc
  dups_anc$TPCA[i] <- tpca
  
  if(tpc > tpca) {p = p+1}
  
}

# have to remove completely identical duplicates bc tpc and maybe even tpca will be 0 
dups_anc <- dups_anc %>%
  filter(TPC > 0 & TPCA > 0)


# extract single copy ortholog pairs that correspond with the dup+anc species pair 
ortho_pairs <- read.csv("./CLOUD/All_Ortholog_Pairs.tsv", sep="")
ortho_pairs <- ortho_pairs %>%
  mutate_at(vars(starts_with("species")), ~ gsub('dana', 'AN', .) %>%
              gsub('dmel', 'ME', .) %>%
              gsub('dmoj', 'MO', .) %>%
              gsub('dper', 'PE', .) %>%
              gsub('dpse', 'PS', .) %>%
              gsub('dvir', 'VI', .) %>%
              gsub('dwil', 'WI', .) %>%
              gsub('dyak', 'YA', .)) %>%
  filter((str_detect(species.x, species_1) & str_detect(species.y, species_2)) |
           (str_detect(species.x, species_2) & str_detect(species.y, species_1)))


#ortho_pairs_rev <- ortho_pairs %>%
#  select(Orthogroup, species.x, YOgn.x, species.y, YOgn.y) %>%
#  mutate(species_pair = paste0(species.y, '_', species.x))
#colnames(ortho_pairs_rev) <- c('Orthogroup','species.y','YOgn.y','species.x','YOgn.x','species_pair')

#ortho_pairs <- rbind(ortho_pairs, ortho_pairs_rev)


# run cloud on each species pair 
func_output <- data.frame()


# get single copy ortholog pairs for the given species pair 
sc_input <- ortho_pairs %>% 
  #filter(species_pair == pair) %>%
  select(YOgn.x,YOgn.y)

# get duplicate+anc species pairs for the given species pair
dup_input <- dups_anc %>%
  #filter(species_pair == pair) %>%
  select(dup_1,TPC,TPCA,dup_2,ancestral_copy)

# merge sc with expression data
exp <- orig_exp
colnames(exp) <- paste0(colnames(orig_exp),'_1')
sc_input <- merge(sc_input, exp, by.x='YOgn.x', by.y='YOgnID_1')

colnames(exp) <- paste0(colnames(orig_exp),'_2')
sc_input <- merge(sc_input, exp, by.x='YOgn.y', by.y='YOgnID_2')

sc_input <- sc_input %>% select(-YOgn.x, -YOgn.y)

write.table(sc_input, paste0('./CLOUD/sc_cloud_input_',pair,'.tsv'))


exp <- orig_exp
orig_dup_input <- dup_input

for (run_1_or_2 in c(1:2)) { # run 1 where dup_1 is parental and run 2 were dup_2 is parental
  dup_input <- orig_dup_input
  
  if(run_1_or_2 == 1){ # merge dups with expression data
    colnames(exp) <- paste0(colnames(orig_exp),'_d1')
    dup_input <- merge(dup_input, exp, by.x='dup_1', by.y='YOgnID_d1')
    
    colnames(exp) <- paste0(colnames(orig_exp),'_d2')
    dup_input <- merge(dup_input, exp, by.x='dup_2', by.y='YOgnID_d2')
  }
  
  if(run_1_or_2 == 2){ # merge dups with expression data but reverse child and parental dups 
    colnames(exp) <- paste0(colnames(orig_exp),'_d2')
    dup_input <- merge(dup_input, exp, by.x='dup_2', by.y='YOgnID_d2')
    
    colnames(exp) <- paste0(colnames(orig_exp),'_d1')
    dup_input <- merge(dup_input, exp, by.x='dup_1', by.y='YOgnID_d1')
  }
  
  # format duplicate gene cloud input
  colnames(exp) <- paste0(colnames(orig_exp),'_anc')
  dup_input <- merge(dup_input, exp, by.x='ancestral_copy', by.y='YOgnID_anc')
  
  dup_input <- dup_input %>%
    select(-dup_1, -dup_2, -ancestral_copy) %>%
    select(TPC, TPCA, everything())
  
  colnames(dup_input) <- c('TPC','TPCA',
                           'eP1','eP2','eP3','eP4','eP5','eP6','eP7','eP8','eP9','eP10','eP11','eP12','ep13','eP14',
                           'eC1','eC2','eC3','eC4','eC5','eC6','eC7','eC8','eC9','eC10','eC11','eC12','eC13','eC14',
                           'eA1','eA2','eA3','eA4','eA5','eA6','eA7','eA8','eA9','eA10','eA11','eA12','eA13','eA14')
  
  
  # log transform dup_input the same way CLOUD does to sc 
  minexp = 1e-4
  errorexp = 1e-5
  m = 14
  for(row_i in 1:nrow(dup_input)) {
    for(col_j in 3:ncol(dup_input)) {
      dup_input[row_i,col_j] <- log10(dup_input[row_i,col_j] + minexp + rnorm(1, 0, errorexp)) # Transform data to log scale, accounting for expression of 0
    }
  }
  
  
  
  write.table(dup_input, paste0('./CLOUD/dup_cloud_input_',pair,'.tsv'))
  
  
  # run cloud
  GenerateTrainingData(m = 14, 
                       Nk = 5, # 50000
                       singlecopy_filename = paste0('./CLOUD/sc_cloud_input_',pair,'.tsv'),
                       duplicate_filename = paste0('./CLOUD/dup_cloud_input_',pair,'.tsv'), 
                       training_prefix = paste0('./CLOUD/First_Try_',pair))

  ClassifierCV(m = 14, 
               batchsize = 5, #50000 
               num_epochs = 5, #500
               log_lambda_min = -5, 
               log_lambda_max = -1, 
               num_lambda = 5, 
               gamma_min = 0, 
               gamma_max = 1, 
               num_gamma = 3,
               training_prefix = paste0('./CLOUD/First_Try_',pair))

  
  read.csv(paste0('./CLOUD/dup_cloud_input_',pair,'.tsv'), sep="") %>%
    mutate(tPC = TPC/TPCA) %>%
    write.table(file = paste0('./CLOUD/dup_cloud_input_2',pair,'.tsv'))
  
  GenerateFeatures(m = 14,
                   singlecopy_filename = paste0('./CLOUD/sc_cloud_input_',pair,'.tsv'),
                   input_filename = paste0('./CLOUD/dup_cloud_input_2',pair,'.tsv'),
                   feature_filename = paste0('./CLOUD/feature_filename_',pair,'.features'))
  
  CLOUDClassify(paste0('./CLOUD/First_Try_',pair),paste0('./CLOUD/feature_filename_',pair))
  
  

  
  # format cloud output
  classes <- read.csv(paste0('./CLOUD/feature_filename_',pair,'.classifications'), sep='')
  
  dup_ids_func <- dups_anc %>%
    #filter(species_pair == pair) %>%
    select(dup_1,dup_2,ancestral_copy) 
  
  if (run_1_or_2 == 1) {
    dups_ids_func_1 <- dup_ids_func %>%
      mutate(func_dup_1_parent = classes$Class)
  }
  
  if (run_1_or_2 == 2){
    dup_ids_func <- dups_ids_func_1 %>% mutate(func_dup_1_child = classes$Class)
    func_output <- rbind(func_output, dup_ids_func)
  }
}

write.table(func_output, file = paste0('./CLOUD/CLOUD_Output_2', pair, '.tsv'))


# model needs to be restarted between runs. or else, dense layers are continuously added





# around 20% of gene trees have
# duplicate copy more similar to an ortholog than its duplicate 
# for example:
#plot(read.tree(text='(YOgnWI05192:2.32249e-06,((YOgnME06467:0.208784,YOgnPS05511:0.282752):1.2146,(YOgnYA01566:0.066968,(YOgnME01797:0.234536,YOgnAN02072:0.099383):0.106142):0.111734):2.32249e-06);'))

#plot(read.tree(text='((YOgnMO06095:0.051567,YOgnVI00812:0.054434):0.0096645,(YOgnAN01058:0.052566,(((YOgnPE05008:0,YOgnPS10129:0):0.036867,(YOgnYA04660:0.002914,YOgnME03065:0.004619):0.027064):0.002351,(YOgnWI07191:0.029138,YOgnYA01355:0.913833):0.000507):0.010438):0.0096645);'))






