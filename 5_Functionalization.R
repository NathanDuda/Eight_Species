


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
dup_anc_exp <- dups %>%  merge(dup_1_exp, ., by.x = 'YOgnID', by.y ='dup_1') %>% rename(YOgnID = 'dup_1')

dup_2_exp <- exp %>% rename_at(-1, ~paste('dup_2_', ., sep = ''))
dup_anc_exp <- dup_anc_exp %>% merge(dup_2_exp, ., by.x = 'YOgnID', by.y ='dup_2') %>% rename(YOgnID = 'dup_2')

anc_exp <- exp %>% rename_at(-1, ~paste('anc_', ., sep = ''))
dup_anc_exp <- dup_anc_exp %>% merge(anc_exp, ., by.x = 'YOgnID', by.y ='ancestral_copy') %>% rename(YOgnID = 'anc')

# get exp for ortholog pairs 
ortho_x_exp <- exp %>% rename_at(-1, ~paste('ortho_x_', ., sep = ''))
ortho_pair_exp <- ortho_pairs %>% merge(ortho_x_exp, ., by.x = 'YOgnID', by.y ='YOgn.x') %>% rename(YOgnID = 'ortho_x')

ortho_y_exp <- exp %>% rename_at(-1, ~paste('ortho_y_', ., sep = ''))
ortho_pair_exp <- ortho_pair_exp %>% merge(ortho_y_exp, ., by.x = 'YOgnID', by.y ='YOgn.y') %>% rename(YOgnID = 'ortho_y') %>%
  select(-species.x, -species.y)

# get combined expression values and after adding, calculate relative expression
dup_1_exp <- raw_exp %>% rename_at(-1, ~paste('dup_1_', ., sep = ''))
dups_combined_exp <- dups %>%  merge(dup_1_exp, ., by.x = 'YOgnID', by.y ='dup_1') %>% rename(YOgnID = 'dup_1')

dup_2_exp <- raw_exp %>% rename_at(-1, ~paste('dup_2_', ., sep = ''))
dups_combined_exp <- dups_combined_exp %>% merge(dup_2_exp, ., by.x = 'YOgnID', by.y ='dup_2') %>% rename(YOgnID = 'dup_2')

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


# use euclidean distance values to classify into functional groups 

strength_coefficient <- 2

func <- ed_values %>%
  rowwise() %>%
  mutate(func = case_when((dup1_a <= cutoff) & (dup2_a <= cutoff) ~ 'conserv',
                          (dup1_a > cutoff) & (dup2_a <= cutoff) ~ 'neo_dup1',
                          (dup1_a <= cutoff & dup2_a > cutoff) ~ 'neo_dup2',
                          (dup1_a > cutoff & dup2_a > cutoff & d1d2_a <= cutoff) ~ 'subfun',
                          (dup1_a > cutoff & dup2_a > cutoff & d1d2_a > cutoff) ~ 'specializ')) %>%
  mutate(func_conservative = 
                case_when((dup1_a <= cutoff*strength_coefficient) & (dup2_a <= cutoff*strength_coefficient) ~ 'conserv',
                          (dup1_a > cutoff*strength_coefficient) & (dup2_a <= cutoff*strength_coefficient) ~ 'neo_dup1',
                          (dup1_a <= cutoff*strength_coefficient & dup2_a > cutoff*strength_coefficient) ~ 'neo_dup2',
                          (dup1_a > cutoff*strength_coefficient & dup2_a > cutoff*strength_coefficient & d1d2_a <= cutoff*strength_coefficient) ~ 'subfun',
                          (dup1_a > cutoff*strength_coefficient & dup2_a > cutoff*strength_coefficient & d1d2_a > cutoff*strength_coefficient) ~ 'specializ')) %>%
  mutate(func_strength =
           case_when(
             func == 'neo_dup1' & func_conservative == 'neo_dup1' ~ 'strong_neo_dup1',
             func == 'neo_dup2' & func_conservative == 'neo_dup2' ~ 'strong_neo_dup2',
             func == 'specializ' & func_conservative == 'specializ' ~ 'strong_specializ',
             func == 'subfun' & func_conservative == 'subfun' ~ 'strong_subfun',
             func == 'conserv' & func_conservative == 'conserv' ~ 'strong_conserv',
           
             func == 'neo_dup1' & func_conservative != 'neo_dup1' ~ 'weak_neo_dup1',
             func == 'neo_dup2' & func_conservative != 'neo_dup2' ~ 'weak_neo_dup2',
             func == 'specializ' & func_conservative != 'specializ' ~ 'weak_specializ',
             func == 'subfun' & func_conservative != 'subfun' ~ 'weak_subfun',
             func == 'conserv' & func_conservative != 'conserv' ~ 'weak_conserv'))

# add dup and anc gene ids to func 
func <- merge(dups,func,by='Orthogroup')

# classify duplicates as having gained a new function if they have expression in a tissue that the ancestral copy was not expressed in
new_exp_neo <- dup_anc_exp
new_exp_neo$func <- NA

for (tissue in tissue_names) {
  new_exp_neo <- new_exp_neo %>%
    mutate(func = case_when(
      select(., paste0('dup_2_',tissue)) != 0 & select(., paste0('dup_1_',tissue)) != 0 & 
        select(., paste0('anc_',tissue)) == 0 ~ 'specializ',
      select(., paste0('dup_1_',tissue)) != 0 & select(., paste0('anc_',tissue)) == 0 ~ 'neo_dup1',
      select(., paste0('dup_2_',tissue)) != 0 & select(., paste0('anc_',tissue)) == 0 ~ 'neo_dup2',
      T ~ func))
}

new_exp_neo <- na.omit(new_exp_neo[c('Orthogroup','func')])
colnames(new_exp_neo) <- c('Orthogroup','new_exp_func')

#matching_rows <- match(func$Orthogroup, new_exp_neo$Orthogroup)
#func$func <- ifelse(!is.na(matching_rows), new_exp_neo$func[matching_rows], func$func)


func <- func %>%
  left_join(., new_exp_neo, by = 'Orthogroup') %>%
  mutate(func = case_when(new_exp_func == 'specializ' ~ 'specializ',
                          func == 'specializ' ~ 'specializ',
                          func == 'neo_dup1' & new_exp_func == 'neo_dup2' ~ 'specializ',
                          func == 'neo_dup2' & new_exp_func == 'neo_dup1' ~ 'specializ',
                          
                          
                          (func == 'conserv' | func == 'subfun') & new_exp_func == 'neo_dup1' ~ 'neo_dup1',
                          (func == 'conserv' | func == 'subfun') & new_exp_func == 'neo_dup2' ~ 'neo_dup2',
                          (func == 'conserv' | func == 'subfun') & new_exp_func == 'specializ' ~ 'specializ',
                          !is.na(func) ~ func
                          
                        ))

### new tissue func counts
how_often_neo_is_newtissue <- func %>%
  mutate(new_tissue = case_when(
                          #func == 'specializ' ~ 'specializ',
                          func == 'specializ' & new_exp_func == 'neo_dup1' ~ 'specializ_newtissue_onedup',
                          func == 'specializ' & new_exp_func == 'neo_dup2' ~ 'specializ_newtissue_onedup',
                          
                          new_exp_func == 'specializ' ~ 'newtissue_specializ',
                          new_exp_func == 'neo_dup1' & new_exp_func == 'neo_dup2' ~ 'newtissue_specializ',
                          
                          #func == 'neo_dup1' ~ 'neo_dup1',
                          #func == 'neo_dup1' ~ 'neo_dup2',
                          func == 'neo_dup2' & new_exp_func == 'neo_dup2' ~ 'newtissue_neo',
                          func == 'neo_dup1' & new_exp_func == 'neo_dup1' ~ 'newtissue_neo'
  ))

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

