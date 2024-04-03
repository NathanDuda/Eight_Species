# Author: Nathan Duda
# Purpose: 
#   This script organizes the output provided by HyPhy and the output provided by MEGA. 


source('./startup.R')

### clean HyPhy output:
# iterate over each hyphy output file to get branch-specific selection values (branches of gene tree of orthogroup)
p = 0 
entire_hyphy_output <- data.frame()

for (file in list.files('./Evolutionary_Rate/HyPhy_absrel_Orthogroup_Output/')){
  
  # import hyphy output json file 
  file <- paste0('./Evolutionary_Rate/HyPhy_absrel_Orthogroup_Output/',file)
  
  # format json file 
  hyphy_output <- tibble(json = list(fromJSON(file))) %>%
    unnest_wider(json, names_sep = '_') %>%
    select(`json_branch attributes`) %>%
    unnest_wider(col = `json_branch attributes`, names_sep = '_') %>%
    unnest_longer(col = `json_branch attributes_0`, names_repair = 'unique') %>%
    unnest_wider(col = `json_branch attributes_0`, names_sep = '_') %>%
    select(`json_branch attributes_0_id`,
           `json_branch attributes_0_Baseline MG94xREV`,
           `json_branch attributes_0_Baseline MG94xREV omega ratio`,
           `json_branch attributes_0_Full adaptive model`,
           `json_branch attributes_0_Full adaptive model (non-synonymous subs/site)`,
           `json_branch attributes_0_Full adaptive model (synonymous subs/site)`) %>%
    filter(str_detect(`json_branch attributes_0_id`, 'YOgn'))
  
  colnames(hyphy_output) <- c('id', 'baseline_MG94xREV', 'baseline_MG94xREV_omega', 'FAM', 'FAM_dn', 'FAM_ds')
  
  entire_hyphy_output <- rbind(entire_hyphy_output, hyphy_output)
  
  # print progress
  p = p + 1
  print(p)
}

# write hyphy output to file 
write.table(entire_hyphy_output, file='./HyPhy_Output.tsv')


# read in one to two orthogroups 
orthogroups <- read.csv2("./Dup_Pair_Orthologs.tsv", sep="")

# format orthogroup dataframe to get orthogroup of each gene 
orthogroup_gn <- orthogroups %>%
  select(-duplicate_pair_species) %>%
  pivot_longer(., cols=2:9, values_to = "Value") %>%
  separate_rows(Value, sep = ", ") %>%
  mutate(name = gsub('_prot', '', name))

# merge to get orthogroups of each gene with their hyphy results 
colnames(orthogroup_gn) <- c('Orthogroup', 'species', 'id')
orthogroups_omega <- merge(orthogroup_gn, entire_hyphy_output, by='id')

# categorize each row as dup or ortho 
orthogroups_omega <- orthogroups_omega %>%
  group_by(Orthogroup, species) %>%
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


