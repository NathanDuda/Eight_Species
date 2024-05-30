



funcs <- read.csv("./CDROM_CLOUD_pseudo_funcs.tsv", sep="")


funcs <- funcs %>%
  mutate(species = case_when(grepl('AN', dup_1) ~ 'dana',
                             grepl('ME', dup_1) ~ 'dmel',
                             grepl('MO', dup_1) ~ 'dmoj',
                             grepl('PE', dup_1) ~ 'dper',
                             grepl('PS', dup_1) ~ 'dpse',
                             grepl('VI', dup_1) ~ 'dvir',
                             grepl('WI', dup_1) ~ 'dwil',
                             grepl('YA', dup_1) ~ 'dyak'))

funcs_counts <- as.data.frame(table(funcs$species, funcs$CLOUD_pseudo))
colnames(funcs_counts) <- c('spec','category','Freq')


ggplot(data = funcs_counts, aes(x="", y=Freq, group=category, fill=category)) +
  geom_bar(width = 1, stat = "identity", position = position_fill()) +
  geom_text(aes(label = Freq), position = position_fill(vjust = 0.5), colour = 'black', size = 3) +
  coord_polar("y", start=0) + 
  facet_grid(.~factor(spec, levels= c('dmel','dyak','dana','dpse','dper','dvir','dmoj','dwil'))) + 
  theme_void() #+
  #scale_fill_manual(values=c("#0079C2", "#FF828B", "808080")) + ###
  #labs(title = 'Mechanism') +
  #theme(legend.position = 'bottom', legend.title = element_blank(), strip.text = element_blank())



