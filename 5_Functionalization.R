


source('./startup.R')


# read in duplicate genes with their ancestral copy 
dups <- read.csv("./Dup_Pairs_Ancestral.tsv", sep="")

# read in ortholog pairs 
ortho_pairs <- read.csv("./Ortholog_Pairs.tsv", sep="")

# read in expression data 
raw_exp <- read.csv("./Expression_Data.tsv", sep="")

# calculate relative expression values 
exp <- raw_exp
rownames(exp) <- exp$YOgnID
exp <- exp %>% select(-YOgnID)
exp <- exp / rowSums(exp)
exp$YOgnID <- rownames(exp)
exp <- exp[, c("YOgnID", setdiff(names(exp), "YOgnID"))]

write.table(exp,'./Relative_Expression.tsv')

# get exp for dups and ancestral
dup_1_exp <- exp %>% rename_at(-1, ~paste('dup_1_', ., sep = ''))
dup_anc_exp <- dups %>%  merge(dup_1_exp, ., by.x = 'YOgnID', by.y ='dup_1') %>% rename(dup_1 = 'YOgnID')

dup_2_exp <- exp %>% rename_at(-1, ~paste('dup_2_', ., sep = ''))
dup_anc_exp <- dup_anc_exp %>% merge(dup_2_exp, ., by.x = 'YOgnID', by.y ='dup_2') %>% rename(dup_2 = 'YOgnID')

anc_exp <- exp %>% rename_at(-1, ~paste('anc_', ., sep = ''))
dup_anc_exp <- dup_anc_exp %>% merge(anc_exp, ., by.x = 'YOgnID', by.y ='ancestral_copy') %>% rename(anc = 'YOgnID')

# get exp for ortholog pairs 
ortho_x_exp <- exp %>% rename_at(-1, ~paste('ortho_x_', ., sep = ''))
ortho_pair_exp <- ortho_pairs %>% merge(ortho_x_exp, ., by.x = 'YOgnID', by.y ='YOgn.x') %>% rename(ortho_x = 'YOgnID')

ortho_y_exp <- exp %>% rename_at(-1, ~paste('ortho_y_', ., sep = ''))
ortho_pair_exp <- ortho_pair_exp %>% merge(ortho_y_exp, ., by.x = 'YOgnID', by.y ='YOgn.y') %>% rename(ortho_y = 'YOgnID') %>%
  select(-species.x, -species.y)

# get combined expression values and after adding, calculate relative expression
dup_1_exp <- raw_exp %>% rename_at(-1, ~paste('dup_1_', ., sep = ''))
dups_combined_exp <- dups %>%  merge(dup_1_exp, ., by.x = 'YOgnID', by.y ='dup_1') %>% rename(dup_1 = 'YOgnID')

dup_2_exp <- raw_exp %>% rename_at(-1, ~paste('dup_2_', ., sep = ''))
dups_combined_exp <- dups_combined_exp %>% merge(dup_2_exp, ., by.x = 'YOgnID', by.y ='dup_2') %>% rename(dup_2 = 'YOgnID')

tissue_names <- c('f_ac','f_dg','f_go','f_hd','f_re','f_tx','f_wb',
                  'm_ac','m_dg','m_go','m_hd','m_re','m_tx','m_wb')
for (tissue in tissue_names) {
  dups_combined_exp <- dups_combined_exp %>%
    rowwise() %>%
    mutate("{tissue}" :=  sum(c_across(ends_with(tissue)), na.rm = TRUE))
}

dups_combined_exp <- dups_combined_exp %>% 
  select(Orthogroup,all_of(tissue_names))

# calculate relative expression values of combined 
dups_combined_exp <- as.data.frame(dups_combined_exp) 
rownames(dups_combined_exp) <- dups_combined_exp$Orthogroup
dups_combined_exp <- dups_combined_exp %>% select(-Orthogroup)
dups_combined_exp <- dups_combined_exp / rowSums(dups_combined_exp)
dups_combined_exp$Orthogroup <- rownames(dups_combined_exp)

# add combined expression values to the dup_anc_exp dataframe  
dups_combined_exp <- dups_combined_exp %>% rename_at(-ncol(dups_combined_exp), ~paste('d1d2_', ., sep = ''))
dup_anc_exp <- merge(dups_combined_exp,dup_anc_exp,by='Orthogroup')

# calculate euclidean distance values 
ed_values <- dup_anc_exp %>%
  select(-dup_1, -dup_2, -anc) %>%
  mutate(dup1_a = sqrt(rowSums((select(., starts_with('dup_1')) - select(., starts_with('anc_'))) ^ 2))) %>%
  mutate(dup2_a = sqrt(rowSums((select(., starts_with('dup_2')) - select(., starts_with('anc_'))) ^ 2))) %>%
  mutate(d1d2_a = sqrt(rowSums((select(., starts_with('d1d2_')) - select(., starts_with('anc_'))) ^ 2))) %>%
  select(Orthogroup, dup1_a, dup2_a, d1d2_a)

sc_ed_values <- ortho_pair_exp %>%
  select(-ortho_x, -ortho_y) %>%
  mutate(sc_ed = sqrt(rowSums((select(., starts_with('ortho_x')) - select(., starts_with('ortho_y'))) ^ 2))) %>%
  select(Orthogroup, sc_ed)
  

# calculate the cutoff ed value 
iqr <- IQR(sc_ed_values$sc_ed) / 2
cutoff <- median(sc_ed_values$sc_ed) + iqr
#cutoff <- cutoff /2


# use euclidean distance values to classify into functional groups 
func <- ed_values %>%
  rowwise() %>%
  mutate(func = case_when((dup1_a <= cutoff) & (dup2_a <= cutoff) ~ 'conserv',
                          (dup1_a > cutoff) & (dup2_a <= cutoff) ~ 'neo_dup1',
                          (dup1_a <= cutoff & dup2_a > cutoff) ~ 'neo_dup2',
                          (dup1_a > cutoff & dup2_a > cutoff & d1d2_a <= cutoff) ~ 'subfun',
                          (dup1_a > cutoff & dup2_a > cutoff & d1d2_a > cutoff) ~ 'specializ'))

# add dup and anc gene ids to func 
func <- merge(dups,func,by='Orthogroup')

# write functionalization results to file
write.table(func, file= 'Dup_Functionalizations.tsv')



