


source('./startup.R')


dups <- read.csv("./Dup_Pairs_Ancestral.tsv", sep="")

orthogroups <- read.csv("./Dup_Pair_Orthologs.tsv", sep="")

exp_rel <- read.csv("./Relative_Expression.tsv", sep="")
exp <- read.csv("./Expression_Data.tsv", sep="")

tau <- read.csv("./Tau.tsv", sep="")

dnds <- read.csv("./orthogroups_omega.tsv", sep="")

longest_transcript <- read.csv("./longest_transcript.tsv", sep="")
seqs <- read.csv("./longest_transcript.tsv", sep="")
seqs <- seqs %>% mutate(gc_content = (str_count(longest_ORF, "G") + str_count(longest_ORF, "C")) / nchar(longest_ORF) * 100)

exon_counts <- read.csv("./Exon_counts.tsv", sep="")

mech <- read.csv("./Dup_Mechanism.tsv", sep="")
mech <- mech[c('Orthogroup','mech')]

func <- read.csv("./Dup_Functionalizations.tsv", sep="")
func <- func[c('Orthogroup','func_strength')]
colnames(func) <- c('Orthogroup','func')

func$func <- gsub('strong_neo_dup1','Strong_Neo',func$func)
func$func <- gsub('strong_neo_dup2','Strong_Neo',func$func)
func$func <- gsub('weak_neo_dup2','Weak_Neo',func$func)
func$func <- gsub('weak_neo_dup1','Weak_Neo',func$func)
func$func <- gsub('weak_specializ','Weak_Spec',func$func)
func$func <- gsub('weak_subfun','Weak_Sub',func$func)
func$func <- gsub('strong_subfun','Strong_Sub',func$func)
func$func <- gsub('strong_conserv','Strong_Cons',func$func)
func$func <- gsub('strong_specializ','Strong_Spec',func$func)





# per orthogroup figures:


# func x mech table 
orthogroup <- merge(dups,func,by='Orthogroup')
orthogroup <- merge(orthogroup,mech,by='Orthogroup')

orthogroup <- orthogroup[orthogroup$Orthogroup %in% dnds$Orthogroup,]

orthogroup$mech <- gsub('RNA_dup_1','RNA',orthogroup$mech)
orthogroup$mech <- gsub('RNA_dup_2','RNA',orthogroup$mech)

func_mech_table <- as.data.frame(unclass(table(orthogroup$mech,orthogroup$func)))
func_mech_table <- func_mech_table[c('Strong_Cons','Weak_Neo','Strong_Neo','Weak_Spec','Strong_Spec','Weak_Sub','Strong_Sub')]

library(knitr)
library(kableExtra)
func_mech_table %>% 
  kable(align = rep('c', 5)) %>%
  column_spec (1:8,border_left = T, border_right = T) %>%
  kable_styling()


# func pie 
func_pie <- orthogroup[c('func')]
func_pie <- as.data.frame(table(func_pie$func))
colnames(func_pie) <- c('func','value')

func_levels <- c('Strong_Cons','Weak_Neo','Strong_Neo','Weak_Spec','Strong_Spec','Weak_Sub','Strong_Sub')
func_colors <- c("#3B7BBD", "#E15F61", "#990000",'#6CBC4D','#267300','#F18244','#cc5200')

#dark_colors  <- c('#004080','#990000','#267300','#cc5200')
#light_colors <- c('#66a3ff','#ff6666','#b3ff66','#ffb366')
#old_colors <- c("#3B7BBD", "#E23F51", "#6CBC4D",'#F18244')

ggplot(func_pie, aes(x="", y=value, group=factor(func,levels=func_levels), fill=factor(func,levels=func_levels))) +
  geom_bar(width = 1, stat = "identity", position = position_fill()) +
  geom_text(aes(label = value), position = position_fill(vjust = 0.5), colour = 'black', size = 3) +
  theme_void() +
  theme(legend.title = element_blank(), panel.background = element_rect(fill = "white", color = "white"), legend.background = element_rect(fill = "white", color = 'white')) +
  coord_polar("y") +
  scale_fill_manual(values=func_colors)

ggsave(filename = './Supplementary_Figures/func_pie.jpg', width = 5, height = 3, device='tiff', dpi=300)




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



dups_ds_func <- as.data.frame(table(dups_ds_func$func,dups_ds_func$Ds_group))

ggplot(dups_ds_func, aes(fill=factor(Var1,levels=func_levels), y=Freq, x=factor(Var2, levels=c('0.125','0.25','0.375','0.5','0.625','0.75','0.875','1','>1')))) + 
  geom_bar(position="fill", stat="identity") +
  scale_fill_manual(values = func_colors)+
  ylab('Percentage') +
  xlab('Ds') +
  theme_classic() +
  theme(legend.title = element_blank())

ggsave(filename = './Supplementary_Figures/ds_grouped_func_barplot.jpg', width = 7, height = 4, device='tiff', dpi=300)



# dup copy figures: 
dups_dnds <- dnds[c('Orthogroup','id','MG94xREV_omega','full_adaptive_model_Dn','full_adaptive_model_Ds','dup_or_ortho')]

# dn, ds, dnds

dups_dnds <- full_join(dups_dnds,func,by='Orthogroup')
dups_dnds <- full_join(dups_dnds,mech,by='Orthogroup')
dups_dnds <- full_join(dups_dnds,dups,by='Orthogroup')
# some one to two orthogroups not matching to dup pair bc not expressed 
# some dup pairs not matching to orthogroups because they dont have at least 4 genes (req by hyphy)



all_dups <- dups_dnds[dups_dnds$dup_or_ortho == 'dup' & !is.na(dups_dnds$dup_or_ortho),]
all_dups$category <- 'All_Dup'
all_ortho <- dups_dnds[dups_dnds$dup_or_ortho == 'ortho' & !is.na(dups_dnds$dup_or_ortho),]
all_ortho$category <- 'All_Ortho'

all_func <- dups_dnds %>%
  filter(!is.na(func) & dup_or_ortho=='dup') %>%
  mutate(category = func)

all_mech <- dups_dnds %>%
  filter(!is.na(mech) & dup_or_ortho=='dup' & mech!='unknown') %>%
  mutate(category = case_when(mech=='RNA_dup_1' & dup_1 == id ~ 'RNA(new)',
                              mech=='RNA_dup_1' & dup_2 == id ~ 'RNA(anc)',
                              mech=='RNA_dup_2' & dup_1 == id ~ 'RNA(anc)',
                              mech=='RNA_dup_2' & dup_2 == id ~ 'RNA(new)',
                              T ~ mech))


all_dnds <- rbind(all_dups,all_ortho,all_func,all_mech)


all_levels <- c('All_Ortho','All_Dup','DNA','RNA(new)','RNA(anc)',func_levels)
all_colors <- c('#778899','#D4E056','#F67280',"#E9AB17","#E9AB17",func_colors)

dnds_plot <- ggplot(all_dnds, aes(x=factor(category, levels = all_levels),y=MG94xREV_omega))+ 
  # facet_grid(cols=vars(variable)) +
  ylab('Dn / Ds') +
  geom_boxplot(color=all_colors) + 
  theme(axis.title.x=element_blank(),panel.background = element_rect(colour = "black", fill=NA, linewidth=1),
        #axis.text.x=element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 20),
        strip.background =element_rect(color='black',fill="white")) +
  ylim(0.000000000000000000000000001,2) +
  #geom_jitter(size=0.001) + 
  stat_compare_means(aes(label = after_stat(p.signif)),method = "t.test", 
                     ref.group = 'All_Ortho',hide.ns = TRUE, size=5) +
  geom_hline(yintercept=median((all_dnds[all_dnds$category=='All_Ortho',])['MG94xREV_omega']$MG94xREV_omega),
             linetype="dashed", color = '#778899',linewidth=0.5) 


dn_plot <- ggplot(all_dnds, aes(x=factor(category, levels = all_levels),y=full_adaptive_model_Dn))+ 
  # facet_grid(cols=vars(variable)) +
  ylab('Dn') +
  geom_boxplot(color=all_colors) + 
  theme(axis.title.x=element_blank(),panel.background = element_rect(colour = "black", fill=NA, size=1),
        #axis.text.x=element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 20),
        strip.background =element_rect(color='black',fill="white")) +
  ylim(0.000000000000000000000000001,2) +
  #geom_jitter(size=0.001) + 
  stat_compare_means(aes(label = after_stat(p.signif)),method = "t.test", 
                     ref.group = 'All_Ortho',hide.ns = TRUE, size=5) +
  geom_hline(yintercept=median((all_dnds[all_dnds$category=='All_Ortho',])['full_adaptive_model_Dn']$full_adaptive_model_Dn),
             linetype="dashed", color = '#778899',linewidth=0.5) 


ds_plot <- ggplot(all_dnds, aes(x=factor(category, levels = all_levels),y=full_adaptive_model_Ds))+ 
  # facet_grid(cols=vars(variable)) +
  ylab('Ds') +
  geom_boxplot(color=all_colors) + 
  theme(axis.title.x=element_blank(),panel.background = element_rect(colour = "black", fill=NA, size=1),
        #axis.text.x=element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 20),
        strip.background =element_rect(color='black',fill="white")) +
  ylim(0.000000000000000000000000001,2)  + 
  stat_compare_means(aes(label = after_stat(p.signif)),method = "t.test", 
                     ref.group = 'All_Ortho',hide.ns = TRUE, size=5) +
  geom_hline(yintercept=median((all_dnds[all_dnds$category=='All_Ortho',])['full_adaptive_model_Ds']$full_adaptive_model_Ds),
             linetype="dashed", color = '#778899',linewidth=0.5) 



ggarrange(dnds_plot,dn_plot,ds_plot,nrow = 3,ncol = 1)

ggsave(filename = './Supplementary_Figures/dnds_dn_ds_boxplot.jpg', width = 11, height = 13)


# tau boxplot
all_tau <- merge(tau,all_dnds,by.x='YOgnID',by.y='id')


ggplot(all_tau, aes(x=factor(category, levels = all_levels),y=tau))+ 
  # facet_grid(cols=vars(variable)) +
  ylab('Tau') +
  geom_boxplot(color = all_colors) + 
  theme(axis.title.x=element_blank(),panel.background = element_rect(colour = "black", fill=NA, size=1),
        #axis.text.x=element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 20),
        strip.background =element_rect(color='black',fill="white")) +
  ylim(0,1) +
  #geom_jitter(size=0.001) + 
  stat_compare_means(aes(label = after_stat(p.signif)),method = "t.test", 
                     ref.group = 'All_Ortho',hide.ns = TRUE, size=5) +
  geom_hline(yintercept=median((all_tau[all_tau$category=='All_Ortho',])['tau']$tau),
             linetype="dashed", color = '#778899',linewidth=0.5) 

ggsave(filename = './Supplementary_Figures/tau_boxplot.jpg', width = 7, height = 4)

write.table(all_tau,'all_plotting_categories.tsv') # REMOVE IF PPI REMOVED OR MOVED INTO THIS SCRIPT

# tissue expression heatmap 

tissue_exp <- pivot_longer(exp_rel, cols = c('f_ac':'m_wb'))

all_tissue_exp <- merge(tissue_exp,all_tau,by='YOgnID')

all_tissue_exp <- all_tissue_exp %>% 
  select(YOgnID, name, value, category) %>% 
  group_by(name,category) %>%
  mutate(value=mean(value)) %>%
  select(name,category,value) %>%
  ungroup() %>%
  distinct()



ggplot(all_tissue_exp, aes(x=factor(category, levels = all_levels), y=name, fill= value)) +
  geom_tile() +
  scale_fill_gradient(high = "dodgerblue4", low = "#9AC0FF",guide = "colorbar") +
  theme(legend.title = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text=element_text(size=8),
        panel.background = element_rect(fill = "white", color = "white"))  

ggsave(filename = './Supplementary_Figures/expression_heatmap.jpg', width = 14, height = 7)


