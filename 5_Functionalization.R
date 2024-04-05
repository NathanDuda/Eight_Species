# Author: Nathan Duda
# Purpose: 
#   This script classifies duplicates into their functional classes using CDROM.

# NOTE:
# from the original CDROM.R script
# remove the chunks of code dealing with plotting densities.
# some of my species pairs only have 1 duplicate pair + ancestral copy occurrence for the species pair
# so error is given when tried to plot density with less than 2 values. 


source('./startup.R')


source('./CDROM/CDROM.R')

dups <- read.csv("./Dup_Pairs_Ancestral.tsv", sep="") %>%
  mutate(dup_species = gsub("[^A-Z]", "", dup_1),
         dup_species = gsub("YO", "", dup_species),
         anc_species = gsub("[^A-Z]", "", ancestral_copy),
         anc_species = gsub("YO", "", anc_species),
         
         species_pair = if_else(dup_species > anc_species,
                                paste0(dup_species, '_', anc_species),
                                paste0(anc_species, '_', dup_species)))


all_cdrom_output <- data.frame()

for (pair in unique(dups$species_pair)) {
  
  species_pair <- strsplit(pair, "_")[[1]]
  species_1 <- species_pair[1]
  species_2 <- species_pair[2]
  
  
  
  # format duplicate genes with their ancestral copy 
  read.csv("./Dup_Pairs_Ancestral.tsv", sep="") %>%
    select(dup_1, dup_2, ancestral_copy) %>%
    filter((grepl(species_1, dup_1) & grepl(species_2, ancestral_copy)) |
             (grepl(species_2, dup_1) & grepl(species_1, ancestral_copy))) %>%
    write.table(., file = './CDROM/dup_file.tsv')
  
  
  # format ortholog pairs 
  orthos <- read.csv("./All_Ortholog_Pairs.tsv", sep="") %>%
    select(YOgn.x, YOgn.y) %>%
    filter((grepl(species_1, YOgn.x) & grepl(species_2, YOgn.y)))
  
  read.csv("./All_Ortholog_Pairs.tsv", sep="") %>%
    select(YOgn.x, YOgn.y) %>%
    filter((grepl(species_1, YOgn.y) & grepl(species_2, YOgn.x))) %>%
    mutate(YOgn.x_old_temp = YOgn.x,
           YOgn.x = YOgn.y,
           YOgn.y = YOgn.x_old_temp) %>%
    select(-YOgn.x_old_temp) %>%
    rbind(., orthos) %>%
    write.table(., file = './CDROM/single_file.tsv')
  
  
  # format expression files
  read.csv("./Expressed_Expression_Data.tsv", sep="") %>%
    filter(grepl(species_1, YOgnID)) %>%
    remove_rownames() %>%
    column_to_rownames('YOgnID') %>%
    write.table(., file = './CDROM/expression_file_1.tsv')
  
  read.csv("./Expressed_Expression_Data.tsv", sep="") %>%
    filter(grepl(species_2, YOgnID)) %>%
    remove_rownames() %>%
    column_to_rownames('YOgnID') %>%
    write.table(., file = './CDROM/expression_file_2.tsv')
  
  # run CDROM 
  CDROM(dupFile = './CDROM/dup_file.tsv',
        singleFile = './CDROM/single_file.tsv',
        exprFile1 = './CDROM/expression_file_1.tsv',
        exprFile2 = './CDROM/expression_file_2.tsv',
        out = './CDROM/Output_',
        PC = FALSE,
        useAbsExpr = FALSE)
  
  
  cdrom_output <- read.delim("./CDROM/Output_1.txt")
  
  all_cdrom_output <- rbind(all_cdrom_output, cdrom_output)

}

write.table(all_cdrom_output, file = './CDROM/All_CDROM_Output.tsv')

all_cdrom_output <- read.csv("C:/Users/17735/Downloads/Eight_Species/CDROM/All_CDROM_Output.tsv", sep="")
func <- all_cdrom_output %>%
  mutate(func = Classification,
         func = case_when(func == 'Conservation' ~ 'cons',
                          func == 'Neofunctionalization(Dup1)' ~ 'neo_dup1',
                          func == 'Neofunctionalization(Dup2)' ~ 'neo_dup2',
                          func == 'Specialization' ~ 'spec',
                          func == 'Subfunctionalization' ~ 'sub')) %>%
  select(-Classification)



## get relative expression levels 
# read in expression data 
raw_exp <- read.csv("./Expressed_Expression_Data.tsv", sep="")

# calculate relative expression values
exp <- raw_exp %>%
  remove_rownames() %>%
  column_to_rownames('YOgnID') %>%
  {. / rowSums(.)} %>%
  rownames_to_column('YOgnID')

# write to file 
write.table(exp,'./Relative_Expression.tsv')


## classify duplicates as having gained a new function if they have expression in
## a tissue that the ancestral copy was not expressed in

dups <- dups %>% select(Orthogroup, dup_1, dup_2, ancestral_copy)

dup_1_exp <- exp %>% rename_at(-1, ~paste('dup_1_', ., sep = ''))
dup_anc_exp <- dups %>%  merge(dup_1_exp, ., by.x = 'YOgnID', by.y ='dup_1') %>% rename(YOgnID = 'dup_1')

dup_2_exp <- exp %>% rename_at(-1, ~paste('dup_2_', ., sep = ''))
dup_anc_exp <- dup_anc_exp %>% merge(dup_2_exp, ., by.x = 'YOgnID', by.y ='dup_2') %>% rename(YOgnID = 'dup_2')

anc_exp <- exp %>% rename_at(-1, ~paste('anc_', ., sep = ''))
dup_anc_exp <- dup_anc_exp %>% merge(anc_exp, ., by.x = 'YOgnID', by.y ='ancestral_copy') %>% rename(YOgnID = 'anc')


new_exp_neo <- dup_anc_exp
new_exp_neo$func <- NA


tissue_names <- c('f_ac','f_dg','f_go','f_hd','f_re','f_tx','f_wb',
                  'm_ac','m_dg','m_go','m_hd','m_re','m_tx','m_wb')


for (tissue in tissue_names) {
  new_exp_neo <- new_exp_neo %>%
    mutate(func = case_when(
      select(., paste0('dup_2_',tissue)) != 0 & select(., paste0('dup_1_',tissue)) != 0 & 
        select(., paste0('anc_',tissue)) == 0 ~ 'spec',
      select(., paste0('dup_1_',tissue)) != 0 & select(., paste0('anc_',tissue)) == 0 ~ 'neo_dup1',
      select(., paste0('dup_2_',tissue)) != 0 & select(., paste0('anc_',tissue)) == 0 ~ 'neo_dup2',
      T ~ func))
}

new_exp_neo <- new_exp_neo %>%
  select(Orthogroup, func) %>%
  na.omit() 

colnames(new_exp_neo) <- c('Orthogroup','new_exp_func')

# attach orthogroups to dups and their functions 
orthogroups <- read.csv("./Dup_Pairs_Ancestral.tsv", sep="") %>%
  select(Orthogroup = Orthogroup, Dup1 = dup_1, Dup2 = dup_2, Ancestor = ancestral_copy, pseudo)
func <- merge(func, orthogroups, by = c('Dup1', 'Dup2', 'Ancestor'))


# combine the original functional classification with my new-tissue functional classification 
func <- func %>%
  left_join(., new_exp_neo, by = 'Orthogroup') %>%
  mutate(func_actual = case_when(
                          new_exp_func == 'spec' ~ 'spec',
                          func == 'spec' ~ 'spec',
                          func == 'neo_dup1' & new_exp_func == 'neo_dup2' ~ 'spec',
                          func == 'neo_dup2' & new_exp_func == 'neo_dup1' ~ 'spec',
                          
                          func == 'cons' & new_exp_func == 'neo_dup1' ~ 'neo_dup1',
                          func == 'cons' & new_exp_func == 'neo_dup2' ~ 'neo_dup2',
                          func == 'cons' & new_exp_func == 'spec' ~ 'spec',
                          
                          func == 'sub' & new_exp_func == 'neo_dup1' ~ 'sub',
                          func == 'sub' & new_exp_func == 'neo_dup2' ~ 'sub',
                          func == 'sub' & new_exp_func == 'spec' ~ 'sub',
                          
                          func == new_exp_func ~ func,
                          is.na(new_exp_func) & !is.na(func) ~ func,
                          is.na(new_exp_func) & is.na(func) ~ pseudo))










new_counts <- as.data.frame((table(how_often_neo_is_newtissue$new_tissue))) %>% 
  pivot_wider(names_from='Var1',values_from='Freq')

func_counts <- as.data.frame(table(how_often_neo_is_newtissue$func)) %>% 
  pivot_wider(names_from='Var1',values_from='Freq')

new_counts$newtissue_specializ / func_counts$specializ # % specialized dups with new tissue expressed in both dups
new_counts$specializ_newtissue_onedup / func_counts$specializ # % specialized dups with new tissue expressed in one dup

new_counts$newtissue_neo / (func_counts$neo_dup1 + func_counts$neo_dup2) # % neo dups with new tissue expressed




###

# write functionalization results to file
write.table(func, file= 'Dup_Functionalizations.tsv')




ggplot(func) +
  #geom_point(aes(x=dup1_a, y=dup2_a, color = func)) +
  geom_point(aes(x=dup1_a, y=dup2_a, color = func_strength)) +
  theme_bw()

