


source('./startup.R')


dups <- read.csv("./Dup_Pairs_Ancestral.tsv", sep="")

tau <- read.csv("./Tau.tsv", sep="")

dnds <- read.csv("./orthogroups_omega.tsv", sep="")

seqs <- read.csv("./longest_transcript.tsv", sep="")
seqs <- seqs %>% mutate(gc_content = (str_count(longest_ORF, "G") + str_count(longest_ORF, "C")) / nchar(longest_ORF) * 100)


mech_exon_counts <- read.csv("./Dup_Mechanism.tsv", sep="")
mech_exon_counts <- mech_exon_counts[c('Orthogroup','dup_1_n_exons','dup_2_n_exons','ancestral_copy_n_exons','mech')]

func <- read.csv("./Dup_Functionalizations.tsv", sep="")
func <- func[c('Orthogroup','func')]
func$func <- gsub('conserv','conserved',func$func)
func$func <- gsub('specializ','specialized',func$func)
func$func <- gsub('subfun','subfunctionalized',func$func)
func$func <- gsub('neo','neofunctionalized',func$func)


# per orthogroup figures:


# func x mech table 
orthogroup <- merge(dups,func,by='Orthogroup')
orthogroup <- merge(orthogroup,mech_exon_counts,by='Orthogroup')

orthogroup <- orthogroup[orthogroup$Orthogroup %in% dnds$Orthogroup,]

func_mech_table <- as.data.frame(unclass(table(orthogroup$mech,orthogroup$func)))
colnames(func_mech_table) <- c('conserved','neo_dup1','neo_dup2','specialized','subfunctionalized')

library(knitr)
library(kableExtra)
func_mech_table %>% 
  kable(align = rep('c', 5)) %>%
  column_spec (1:6,border_left = T, border_right = T) %>%
  kable_styling()


# func pie 
func_pie <- orthogroup[c('func')]

func_pie$func <- gsub('neofunctionalized_dup1','neofunctionalized',func_pie$func)
func_pie$func <- gsub('neofunctionalized_dup2','neofunctionalized',func_pie$func)

func_pie <- as.data.frame(table(func_pie$func))

colnames(func_pie) <- c('func','value')

ggplot(func_pie, aes(x="", y=value, group=func, fill=func)) +
  geom_bar(width = 1, stat = "identity", position = position_fill()) +
  geom_text(aes(label = value), position = position_fill(vjust = 0.5), colour = 'black', size = 3) +
  theme_void() +
  theme(legend.title = element_blank(), panel.background = element_rect(fill = "white", color = "white"), legend.background = element_rect(fill = "white", color = 'white')) +
  coord_polar("y") +
  scale_fill_manual(values=c("#3B7BBD", "#E23F51", "#6CBC4D",'#F18244'))

ggsave(filename = './Plots/func_pie.jpg', width = 5, height = 3, device='tiff', dpi=300)


# mech pie 
mech_pie <- orthogroup[c('mech')]

mech_pie$mech <- gsub('_dup_1','',mech_pie$mech)
mech_pie$mech <- gsub('_dup_2','',mech_pie$mech)

mech_pie <- as.data.frame(table(mech_pie$mech))
colnames(mech_pie) <- c('mech','value')

ggplot(mech_pie, aes(x="", y=value, group=mech, fill=mech)) +
  geom_bar(width = 1, stat = "identity", position = position_fill()) +
  geom_text(aes(label = value), position = position_fill(vjust = 0.5), colour = 'black', size = 3) +
  theme_void() +
  theme(legend.title = element_blank(), panel.background = element_rect(fill = "white", color = "white"), legend.background = element_rect(fill = "white", color = 'white')) +
  coord_polar("y") +
  scale_fill_manual(values=c("#F67280", "#E9AB17", "#1E90FF"))

ggsave(filename = './Plots/mech_pie.jpg', width = 5, height = 3, device='tiff', dpi=300)


# ds group func barplot 
dups_ds_func <- dnds %>%
  mutate(Ds_group=case_when(full_adaptive_model_Ds <= 0.125 & full_adaptive_model_Ds > 0 ~ '0.125',
                            full_adaptive_model_Ds <= 0.25 & full_adaptive_model_Ds > 0.125 ~ '0.25',
                            full_adaptive_model_Ds <= 0.375 & full_adaptive_model_Ds > 0.25 ~ '0.375',
                            full_adaptive_model_Ds <= 0.5 & full_adaptive_model_Ds > 0.375 ~ '0.5',
                            full_adaptive_model_Ds <= 0.625 & full_adaptive_model_Ds > 0.5 ~ '0.625',
                            full_adaptive_model_Ds <= 0.75 & full_adaptive_model_Ds > 0.625 ~ '0.75',
                            full_adaptive_model_Ds <= 0.875 & full_adaptive_model_Ds > 0.75 ~ '0.875',
                            full_adaptive_model_Ds <= 1 & full_adaptive_model_Ds > 0.875 ~ '1',
                            full_adaptive_model_Ds > 1 ~ '>1',))

dups_ds_func <- merge(func,dups_ds_func,by='Orthogroup')

dups_ds_func$func <- gsub('neofunctionalized_dup1','neofunctionalized',dups_ds_func$func)
dups_ds_func$func <- gsub('neofunctionalized_dup2','neofunctionalized',dups_ds_func$func)

dups_ds_func <- as.data.frame(table(dups_ds_func$func,dups_ds_func$Ds_group))

ggplot(dups_ds_func, aes(fill=Var1, y=Freq, x=factor(Var2, levels=c('0.125','0.25','0.375','0.5','0.625','0.75','0.875','1','>1')))) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = c("#3B7BBD",'#E23F51','#6CBC4D','#F18244'))+
  ylab('Percentage') +
  xlab('Ds') +
  theme_classic() +
  theme(legend.title = element_blank())

ggsave(filename = './Plots/ds_grouped_func_barplot.jpg', width = 7, height = 4, device='tiff', dpi=300)



# dup copy figures: 
dups_dnds <- dnds[c('Orthogroup','id','MG94xREV_omega','full_adaptive_model_Dn','full_adaptive_model_Ds','dup_or_ortho')]

# dn, ds, dnds

dups_dnds <- full_join(dups_dnds,func,by='Orthogroup')
dups_dnds <- full_join(dups_dnds,mech_exon_counts,by='Orthogroup')
dups_dnds <- full_join(dups_dnds,dups,by='Orthogroup')
# some one to two orthogroups not matching to dup pair bc not expressed 
# some dup pairs not matching to orthogroups because they dont have at least 4 genes (req by hyphy)











########
# OLD SCRIPT:


alldup_exon_plot$func_copy <- 'duplicate_pairs'

exon_plot <- rbind(exon_plot,alldup_exon_plot)

# add boxplot for all orthologs
ortho_exon_plot <- orthos[,c(5,12)]

ortho_exon_plot$func_copy <- 'ortholog_pairs'
ortho_exon_plot$mech_copy <- 'ortholog_pairs'

exon_plot <- rbind(exon_plot,ortho_exon_plot)

exon_plot <- melt(exon_plot,id.vars = c('mech_copy','func_copy'))
exon_plot <- exon_plot %>% 
  mutate(func_copy = case_when(variable == 'exon_diff' & func_copy == 'not_neo_copy' ~ 'neofunctionalized',
                               variable == 'exon_diff' & func_copy == 'neo_copy' ~ 'neofunctionalized',
                               TRUE ~ func_copy))



# add mechanism boxes
new_RNA <- exon_plot[exon_plot$mech_copy=='new_RNA',]
new_RNA$func_copy <- 'new_RNA'

old_RNA <- exon_plot[exon_plot$mech_copy=='old_RNA',]
old_RNA$func_copy <- 'old_RNA'

DNA <- exon_plot[exon_plot$mech_copy=='DNA',]
DNA$func_copy <- 'DNA'

exon_plot <- rbind(exon_plot,new_RNA,old_RNA,DNA)
#exon_plot <- exon_plot[,-c(1)]
#


exon_plot$variable <- gsub('exons','Combined',exon_plot$variable)
exon_plot$variable <- gsub('exon_diff','Difference',exon_plot$variable)

exon_plot_diff <- exon_plot[exon_plot$variable == 'Difference',]
exon_plot <- exon_plot[exon_plot$variable == 'Combined',]

exon_plot_diff$variable <- 'Difference'
exon_plot$variable <- 'Combined'



func_levels <- c('ortholog_pairs','duplicate_pairs','DNA','new_RNA','old_RNA','conserved','neo_copy','not_neo_copy','specialized','subfunctionalized')
exon_plot <- ggplot(exon_plot, aes(x=factor(func_copy, levels = func_levels),y=value))+ 
  # facet_grid(cols=vars(variable)) +
  ylab('Exon Amount') +
  geom_boxplot(color=c('#778899','#D4E056','#F67280',"#E9AB17","#E9AB17","#3B7BBD",'#E23F51','#E23F51','#6CBC4D','#F18244')) + 
  theme(axis.title.x=element_blank(),panel.background = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 20),
        strip.background =element_rect(color='black',fill="white")) + 
  stat_compare_means(aes(label = after_stat(p.signif)),method = "t.test", 
                     ref.group = 'ortholog_pairs',hide.ns = TRUE,size=6)  +
  geom_hline(yintercept=median((exon_plot[exon_plot$func_copy=='ortholog_pairs',])[4]$value),
             linetype="dashed", color = '#778899',linewidth=0.5)



exon_plot_diff$func_copy <- gsub('not_neo_copy','neofunctionalized',exon_plot_diff$func_copy)
exon_plot_diff$func_copy <- gsub('neo_copy','neofunctionalized',exon_plot_diff$func_copy)

exon_plot_diff$func_copy <- gsub('new_RNA','RNA',exon_plot_diff$func_copy)
exon_plot_diff$func_copy <- gsub('old_RNA','RNA',exon_plot_diff$func_copy)


diff_labels <- c('Ortho','All_Dup','DNA','RNA','Cons','Neo','Spec','Sub')

func_levels_diff <- c('ortholog_pairs','duplicate_pairs','DNA','RNA','conserved','neofunctionalized','specialized','subfunctionalized')
exon_plot_diff <- ggplot(exon_plot_diff, aes(x=factor(func_copy, levels = func_levels_diff),y=value))+ 
  # facet_grid(cols=vars(variable)) +
  geom_boxplot(color=c('#778899','#D4E056','#F67280',"#E9AB17","#3B7BBD",'#E23F51','#6CBC4D','#F18244')) +  
  theme(axis.title.x=element_blank(),panel.background = element_rect(colour = "black", fill=NA, size=1),
        #    axis.text.x=element_blank(),
        axis.text.x = element_text(size = 9),
        axis.text.y=element_blank(),
        axis.title.y=element_blank(),strip.background =element_rect(color='black',fill="white")) + 
  stat_compare_means(aes(label = after_stat(p.signif)),method = "t.test", 
                     ref.group = 'ortholog_pairs',hide.ns = TRUE,size=6)+
  ylab('')  +
  geom_hline(yintercept=median((exon_plot_diff[exon_plot_diff$func_copy=='ortholog_pairs',])[4]$value),
             linetype="dashed", color = '#778899',linewidth=0.5) +
  scale_x_discrete(labels=diff_labels)








# N ORTHOLOGS AND DS AND FUNC 













