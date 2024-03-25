


source('./startup.R')

source("./CLOUD/CLOUD.r")




library(tensorflow)
library(keras) 
library(mvtnorm)


orthogroups <- read.delim2("./OrthoFinder_Output/Results_Jan01/Orthogroups/Orthogroups.tsv", na.strings = '')

og_with_tree <- orthogroups %>%
  mutate(n_na = rowSums(is.na(across(everything())))) %>%
  filter(n_na <= 5)


dups_anc <- read.csv("./Dup_Pairs_Ancestral.tsv", sep="")
dups_anc <- dups_anc %>%
  filter(Orthogroup %in% og_with_tree$Orthogroup)



#t <- dups_anc %>%
#  mutate(across(everything(), ~gsub("\\d|YOgn", "", .))) %>%
#  select(-Orthogroup)
#table(t$dup_1,t$ancestral_copy)
#d <- as.data.frame(table(t$dup_1,t$ancestral_copy))


dups_anc <- dups_anc %>%
  mutate(dup_species = gsub("\\d|YOgn", "", dup_1),
         anc_species = gsub("\\d|YOgn", "", ancestral_copy)) %>%
  filter(case_when(dup_species == 'YA' & anc_species == 'ME' ~ T,
                   dup_species == 'WI' & anc_species == 'PS' ~ T,
                   dup_species == 'PS' & anc_species == 'PE' ~ T,
                   dup_species == 'PE' & anc_species == 'PS' ~ T,
                   dup_species == 'MO' & anc_species == 'VI' ~ T,
                   dup_species == 'AN' & anc_species == 'YA' ~ T,
                   dup_species == 'VI' & anc_species == 'MO' ~ T,
                   dup_species == 'ME' & anc_species == 'YA' ~ T)) %>%
  mutate(species_pair = if_else(dup_species > anc_species,
                                paste0(dup_species,'_',anc_species),
                                paste0(anc_species,'_',dup_species)))




dups_anc$TPC <- NA
dups_anc$TPCA <- NA

p = 0
for(i in 1:nrow(dups_anc)) {
  row <- dups_anc[i,]
  
  OG <- row$Orthogroup
  dup1 <- row$dup_1
  dup2 <- row$dup_2
  anc <- row$ancestral_copy
  
  
  tree <- read.tree(paste0('./OrthoFinder_Output/Results_Jan01/Gene_Trees/',OG,'_tree.txt'))
  
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


ortho_pairs <- read.csv("./All_Ortholog_Pairs.tsv", sep="")
ortho_pairs <- ortho_pairs %>%
  mutate_at(vars(starts_with("species")), ~ gsub('dana', 'AN', .) %>%
                                            gsub('dmel', 'ME', .) %>%
                                            gsub('dmoj', 'MO', .) %>%
                                            gsub('dper', 'PE', .) %>%
                                            gsub('dpse', 'PS', .) %>%
                                            gsub('dvir', 'VI', .) %>%
                                            gsub('dwil', 'WI', .) %>%
                                            gsub('dyak', 'YA', .)) %>%
  mutate(species_pair = if_else(species.x > species.y,
                                paste0(species.x, '_', species.y),
                                paste0(species.y, '_', species.x)))


orig_exp <- read.csv("./Expression_Data.tsv", sep="")

t <- orig_exp %>% # EXPRESSION SHOULD BE LOG TRANSFORMEDDDDDDDDDDDDDDDDDDDDDDDDDDDDDDD
  remove_rownames() %>%
  column_to_rownames('YOgnID') %>%
  log10() %>%
  rownames_to_column()


for (pair in unique(dups_anc$species_pair)) {
  print(pair)
  
  sc_input <- ortho_pairs %>% 
    filter(species_pair == pair) %>%
    select(YOgn.x,YOgn.y)
  
  dup_input <- dups_anc %>%
    filter(species_pair == pair) %>%
    select(dup_1,TPC,TPCA,dup_2,ancestral_copy)
  
  # merge sc with expression data
  exp <- orig_exp
  colnames(exp) <- paste0(colnames(orig_exp),'_1')
  sc_input <- merge(sc_input, exp, by.x='YOgn.x', by.y='YOgnID_1')
  
  colnames(exp) <- paste0(colnames(orig_exp),'_2')
  sc_input <- merge(sc_input, exp, by.x='YOgn.y', by.y='YOgnID_2')
  
  sc_input <- sc_input %>% select(-YOgn.x, -YOgn.y)
  
  write.table(sc_input, './CLOUD/sc_cloud_input.tsv')
  
  
  # merge dups with expression data
  
  exp <- orig_exp
  colnames(exp) <- paste0(colnames(orig_exp),'_d1')
  dup_input <- merge(dup_input, exp, by.x='dup_1', by.y='YOgnID_d1')
  
  colnames(exp) <- paste0(colnames(orig_exp),'_d2')
  dup_input <- merge(dup_input, exp, by.x='dup_2', by.y='YOgnID_d2')
  
  colnames(exp) <- paste0(colnames(orig_exp),'_anc')
  dup_input <- merge(dup_input, exp, by.x='ancestral_copy', by.y='YOgnID_anc')
  
  dup_input <- dup_input %>%
    select(-dup_1, -dup_2, -ancestral_copy) %>%
    select(TPC, TPCA, everything())
  
  write.table(dup_input, './CLOUD/dup_cloud_input.tsv')
  
  
  GenerateTrainingData(m = 14, 
                       Nk = 50,  # change to 50000
                       singlecopy_filename = './CLOUD/sc_cloud_input.tsv', 
                       duplicate_filename = './CLOUD/dup_cloud_input.tsv', 
                       training_prefix = './CLOUD/First_Try_training')
  
  
  ClassifierCV(m = 14, 
               batchsize = 50, # change to 50000
               num_epochs = 50, # change to 500
               log_lambda_min = -5, 
               log_lambda_max = -1, 
               num_lambda = 5, 
               gamma_min = 0, 
               gamma_max = 1, 
               num_gamma = 3, 
               training_prefix = './CLOUD/First_Try_training')
  
  
  PredictorCV(7, 50, 50, -5, -1, 5, 0, 1, 3, './CLOUD/First_Try_training')
  
  
  
  colnames(dup_input) <- c('tPC','tPCA',
                           'eP1','eP2','eP3','eP4','eP5','eP6','eP7','eP8','eP9','eP10','eP11','eP12','ep13','eP14',
                           'eC1','eC2','eC3','eC4','eC5','eC6','eC7','eC8','eC9','eC10','eC11','eC12','eC13','eC14',
                           'eA1','eA2','eA3','eA4','eA5','eA6','eA7','eA8','eA9','eA10','eA11','eA12','eA13','eA14')
  write.table(dup_input, './CLOUD/dup_cloud_input.tsv')
  
  
  #GenerateFeatures(m = 14,
  #                 singlecopy_filename = './CLOUD/sc_cloud_input.tsv',
  #                 input_filename = './CLOUD/dup_cloud_input.tsv',
  #                 feature_filename = './CLOUD/First_Try_testing')
  
  
  GenerateTrainingData(14, 
                       100, # CHANGE THIS NUMBERRRRRRRRRRRRRRRRRRRRRRRRRRR
                       './CLOUD/sc_cloud_input.tsv',
                       'CLOUD/dup_cloud_input.tsv', 
                       './CLOUD/First_Try_testing')
  
  
  
  CLOUDClassify(training_prefix = './CLOUD/First_Try_training', 
                testing_prefix = './CLOUD/First_Try_testing')
  
  
  
  
}


# around 20% of gene trees have
# duplicate copy more similar to an ortholog than its duplicate 
# for example:
plot(read.tree(text='(YOgnWI05192:2.32249e-06,((YOgnME06467:0.208784,YOgnPS05511:0.282752):1.2146,(YOgnYA01566:0.066968,(YOgnME01797:0.234536,YOgnAN02072:0.099383):0.106142):0.111734):2.32249e-06);'))


plot(read.tree(text='((YOgnMO06095:0.051567,YOgnVI00812:0.054434):0.0096645,(YOgnAN01058:0.052566,(((YOgnPE05008:0,YOgnPS10129:0):0.036867,(YOgnYA04660:0.002914,YOgnME03065:0.004619):0.027064):0.002351,(YOgnWI07191:0.029138,YOgnYA01355:0.913833):0.000507):0.010438):0.0096645);'))



ortho_pairs <- read.csv("./All_Ortholog_Pairs.tsv", sep="")













