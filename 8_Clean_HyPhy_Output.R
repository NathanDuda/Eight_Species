
source('./startup.R')


# iterate over each hyphy output file to get 
n = 0 
entire_hyphy_output <- data.frame()

for (file in list.files('./Evolutionary_Rate/HyPhy_busted_Orthogroup_Output/')){
  
  file <- paste0('./Evolutionary_Rate/HyPhy_busted_Orthogroup_Output/',file)
  
  hyphy_output <- jsonlite::fromJSON(file, flatten = TRUE)
  hyphy_output <- hyphy_output['branch attributes']
  
  hyphy_output <- hyphy_output %>%
    as.data.frame() %>%
    t() %>%
    as.data.frame() %>%
    mutate(id = rownames(.)) %>%
    filter(!(grepl('by.site',id) | grepl('Posterior.prob',id) | grepl('original.name',id))) %>%
    filter(grepl('YOgn',id)) %>%
    mutate(id = gsub('branch.attributes.0.','',id)) %>%
    mutate(id = gsub('.with.separate.rates.for.branch.sets','',id)) %>%
    mutate(id = gsub('Nucleotide.GTR','Nucleotide_GTR',id)) %>%
    separate(id, sep='\\.', into = c('id','attribute')) %>%
    rename_with(~ 'value', 1)  %>%
    select(id,attribute,value) %>%
    pivot_wider(names_from='attribute', values_from ='value')
  
  
  # add contrained column if missing. Missing when:
  # no evidence for selection found in unconstrained, so no need to test for selection with constrained 
  if (!'constrained' %in% names(hyphy_output)) {hyphy_output$constrained <- NA}
  
  entire_hyphy_output <- rbind(entire_hyphy_output,hyphy_output)
  
  n = n + 1
  print(n)
}

write.table(entire_hyphy_output,file='./HyPhy_Output.tsv')


orthogroups <- read.csv2("./Dup_Pair_Orthologs.tsv", sep="")

orthogroup_gn <- orthogroups %>%
  select(-duplicate_pair_species) %>%
  pivot_longer(., cols=2:9, values_to = "Value") %>%
  separate_rows(Value, sep = ", ") %>%
  mutate(name = gsub('_prot','',name))

colnames(orthogroup_gn) <- c('Orthogroup','species','id')
orthogroups_omega <- merge(orthogroup_gn,entire_hyphy_output,by='id')


write.table(orthogroups_omega, file = './orthogroups_omega.tsv')


#



t <- orthogroups_omega %>%
  group_by(Orthogroup,species) %>%
  mutate(dup_or_ortho = case_when(n() ==2 ~ 'dup',
                                  n()==1 ~ 'ortho'))



t$MG94xREV <- as.numeric(t$MG94xREV)
t$unconstrained <- as.numeric(t$unconstrained)

ggplot(t,aes(x=dup_or_ortho,y=unconstrained)) +
  geom_boxplot() +
  #scale_y_log10() +
  ylim(0.0000000000001,4)
  theme_bw()








