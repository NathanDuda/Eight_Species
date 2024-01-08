


source('./startup.R')


### clean HyPhy output:

# iterate over each hyphy output file to get branch-specific selection values (branches of gene tree of orthogroup)
p = 0 
entire_hyphy_output <- data.frame()

for (file in list.files('./Evolutionary_Rate/HyPhy_absrel_Orthogroup_Output/')){
  
  # import hyphy output json file 
  file <- paste0('./Evolutionary_Rate/HyPhy_absrel_Orthogroup_Output/',file)
  hyphy_output <- jsonlite::fromJSON(file, flatten = TRUE)
  hyphy_output <- hyphy_output[["branch attributes"]][['0']]
  
  # remove the lists that are only present for some nodes 
  hyphy_output <- hyphy_output[!names(hyphy_output) %in% c("posterior", "original name")]
  hyphy_output <- map(hyphy_output, ~ .x[!names(.x) %in% c("posterior", "original name")])
  
  # duplicate rows to 'Rate Distributions' so that all have four rows  
  for (n in 1:length(hyphy_output)) {
    if (nrow(hyphy_output[[n]][['Rate Distributions']]) == 1){  
      row_to_duplicate <- hyphy_output[[n]][['Rate Distributions']][1, ]
      for (i in c(1,2,3)) {hyphy_output[[n]][['Rate Distributions']] <- rbind(hyphy_output[[n]][['Rate Distributions']], row_to_duplicate)}}
    if (nrow(hyphy_output[[n]][['Rate Distributions']]) == 2){  
      row_to_duplicate <- hyphy_output[[n]][['Rate Distributions']][1, ]
      for (i in c(1,2)) {hyphy_output[[n]][['Rate Distributions']] <- rbind(hyphy_output[[n]][['Rate Distributions']], row_to_duplicate)}}  
    if (nrow(hyphy_output[[n]][['Rate Distributions']]) == 3){  
      row_to_duplicate <- hyphy_output[[n]][['Rate Distributions']][1, ]
      hyphy_output[[n]][['Rate Distributions']] <- rbind(hyphy_output[[n]][['Rate Distributions']], row_to_duplicate)}  
  }
  
  # format hyphy output 
  hyphy_output <- hyphy_output %>%
    
    # make lists into data frame 
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    
    # get names of ids and attributes 
    mutate(id = rownames(.)) %>%
    filter(!(grepl('original.name',id) | grepl('posterior',id))) %>%
    filter(grepl('YOgn',id)) %>%
    mutate(id = gsub('X0.','',id)) %>%
    mutate(id = str_replace(id, "\\.", "_")) %>%
    separate(id, sep='_', into = c('id','attribute')) %>%
    rename_all(~ c('value_1','value_2','value_3','value_4','id','attribute')) %>%
    
    # only keep the extra values in 'Rate.Distributions' because all others were just duplicated
    mutate(value_1 = as.numeric(value_1)) %>%
    mutate(value_2 = case_when(attribute %in% c('Rate.Distributions.1', 'Rate.Distributions.2') ~ as.numeric(value_2), T ~ NA)) %>%
    mutate(value_3 = case_when(attribute %in% c('Rate.Distributions.1', 'Rate.Distributions.2') ~ as.numeric(value_3), T ~ NA)) %>%
    mutate(value_4 = case_when(attribute %in% c('Rate.Distributions.1', 'Rate.Distributions.2') ~ as.numeric(value_4), T ~ NA)) %>%
    
    # pivot wider to get the values in separate columns for each gene which is now a row 
    pivot_wider(names_from='attribute', values_from = c('value_1','value_2','value_3','value_4')) %>%
    
    # remove entire columns that are NA. Those that were duplicated before 
    select(which(!apply(is.na(.), 2, all)))
  
  # add the output part to the entire dataframe 
  entire_hyphy_output <- rbind(entire_hyphy_output,hyphy_output)
  
  # print progress
  p = p + 1
  print(p)
}

# write hyphy output to file 
write.table(entire_hyphy_output,file='./HyPhy_Output.tsv')


# read in one to two orthogroups 
orthogroups <- read.csv2("./Dup_Pair_Orthologs.tsv", sep="")

# format orthogroup dataframe to get orthogroup of each gene 
orthogroup_gn <- orthogroups %>%
  select(-duplicate_pair_species) %>%
  pivot_longer(., cols=2:9, values_to = "Value") %>%
  separate_rows(Value, sep = ", ") %>%
  mutate(name = gsub('_prot','',name))

# merge to get orthogroups of each gene with their hyphy results 
colnames(orthogroup_gn) <- c('Orthogroup','species','id')
orthogroups_omega <- merge(orthogroup_gn,entire_hyphy_output,by='id')

colnames(orthogroups_omega) <- 
  c('id','Orthogroup','species','MG94xREV','MG94xREV_omega','corrected_pval',
  'full_adaptive_model','full_adaptive_model_Dn','full_adaptive_model_Ds','LRT','Nuc_GTR','rate_distributions_1',
  'rate_distributions_2','rate_classes','uncorrected_pval','rate_distributions_3','rate_distributions_4',
  'rate_distributions_5','rate_distributions_6','rate_distributions_7','rate_distributions_8')

# categorize each row as dup or ortho 
orthogroups_omega <- orthogroups_omega %>%
  group_by(Orthogroup,species) %>%
  mutate(dup_or_ortho = case_when(n() == 2 ~ 'dup',
                                  n() == 1 ~ 'ortho'))

# write hyphy results for each orthogroup to file 
write.table(orthogroups_omega, file = './orthogroups_omega.tsv')



### clean MEGA output:
library(phytools)
library(ape)

# read in the newick trees
trees <- ape::read.tree('./MEGA_Phylogeny/MEGA_consensus_Newicks.txt')
trees <- ape::read.tree('./MEGA_Phylogeny/MEGA_Newicks.txt')

# make sure all trees have 8 tips. Nothing should print
for (tree in trees){if (length(tree$tip.label) != 8) {print(tree['tip.label'])}}

# calculate the least squares consensus tree 
consensus_tree <- ls.consensus(trees)

# add root by using midpoint as root
consensus_tree <- midpoint.root(consensus_tree)

# rename the tips 
consensus_tree[["tip.label"]] <- c('dyak','dmel','dmoj','dvir','dper','dpse','dwil','dana')

# plot the phylogeny
plotTree(consensus_tree)

# save the phylogeny to a file 
jpeg("./Plots/Phylogeny.jpg", height = 500, width =850)
plotTree(consensus_tree)
dev.off()

write.tree(consensus_tree, file = "./MEGA_Tree.nwk")


