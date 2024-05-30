



funcs <- read.csv("./CDROM_CLOUD_pseudo_funcs.tsv", sep="")

ds <- read.csv("./orthogroups_omega.tsv", sep="") 

ds <- ds %>%
  mutate(ds_group = case_when(FAM_ds <= 0.125 & FAM_ds > 0 ~ '0.125',
                              FAM_ds <= 0.25 & FAM_ds > 0.125 ~ '0.25',
                              FAM_ds <= 0.375 & FAM_ds > 0.25 ~ '0.375',
                              FAM_ds <= 0.5 & FAM_ds > 0.375 ~ '0.5',
                              FAM_ds <= 0.625 & FAM_ds > 0.5 ~ '0.625',
                              FAM_ds <= 0.75 & FAM_ds > 0.625 ~ '0.75',
                              FAM_ds <= 0.875 & FAM_ds > 0.75 ~ '0.875',
                              FAM_ds <= 1 & FAM_ds > 0.875 ~ '1',
                              FAM_ds > 1 ~ '>1',))

ds <- ds %>%
  select(Orthogroup, id, FAM_ds, ds_group)

t <- funcs %>%
  select(-CDROM_func, -CLOUD_func, -pseudo) %>%
  pivot_longer(values_to = 'id', cols = c(dup_1, dup_2, anc)) %>%
  left_join(., ds, by = 'id') 


# ds group func barplot 



dups_ds_func <- as.data.frame(table(t$CLOUD_pseudo, t$ds_group))

ggplot(dups_ds_func, aes(fill=Var1, y=Freq, x=factor(Var2, levels=c('0.125','0.25','0.375','0.5','0.625','0.75','0.875','1','>1')))) + 
  geom_bar(position="fill", stat="identity") +
  #scale_fill_manual(values = c("#3B7BBD",'#E23F51','#6CBC4D','#F18244'))+
  ylab('Percentage') +
  xlab('Ds') +
  theme_classic() +
  theme(legend.title = element_blank())

ggsave(filename = './Plots/ds_grouped_func_barplot.jpg', width = 7, height = 4, device='tiff', dpi=300)



