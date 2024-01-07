


source('./startup.R')


dups <- read.csv("./Dup_Pairs_Ancestral.tsv", sep="")

exp <- read.csv("./Relative_Expression.tsv", sep="")

tau <- read.csv("./Tau.tsv", sep="")

dnds <- read.csv("./orthogroups_omega.tsv", sep="")

seqs <- read.csv("./longest_transcript.tsv", sep="")
seqs <- seqs %>% mutate(gc_content = (str_count(longest_ORF, "G") + str_count(longest_ORF, "C")) / nchar(longest_ORF) * 100)

exon_counts <- read.csv("./Exon_counts.tsv", sep="")

mech <- read.csv("./Dup_Mechanism.tsv", sep="")
mech <- mech[c('Orthogroup','mech')]

func <- read.csv("./Dup_Functionalizations.tsv", sep="")
func <- func[c('Orthogroup','func')]
func$func <- gsub('conserv','conserved',func$func)
func$func <- gsub('specializ','specialized',func$func)
func$func <- gsub('subfun','subfunctionalized',func$func)
func$func <- gsub('neo','neofunctionalized',func$func)


# per orthogroup figures:


# func x mech table 
orthogroup <- merge(dups,func,by='Orthogroup')
orthogroup <- merge(orthogroup,mech,by='Orthogroup')

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
  mutate(category = case_when(func=='neofunctionalized_dup1' & dup_1 == id ~ 'Neo(new)',
                              func=='neofunctionalized_dup1' & dup_2 == id ~ 'Neo(anc)',
                              func=='neofunctionalized_dup2' & dup_1 == id ~ 'Neo(anc)',
                              func=='neofunctionalized_dup2' & dup_2 == id ~ 'Neo(new)',
                              func=='conserved' ~ 'Cons',
                              func=='specialized' ~ 'Spec',
                              func=='subfunctionalized' ~ 'Sub'))

all_mech <- dups_dnds %>%
  filter(!is.na(mech) & dup_or_ortho=='dup' & mech!='unknown') %>%
  mutate(category = case_when(mech=='RNA_dup_1' & dup_1 == id ~ 'RNA(new)',
                              mech=='RNA_dup_1' & dup_2 == id ~ 'RNA(anc)',
                              mech=='RNA_dup_2' & dup_1 == id ~ 'RNA(anc)',
                              mech=='RNA_dup_2' & dup_2 == id ~ 'RNA(new)',
                              T ~ mech))


all_dnds <- rbind(all_dups,all_ortho,all_func,all_mech)


func_levels <- c('All_Ortho','All_Dup','DNA','RNA(new)','RNA(anc)','Cons','Neo(new)','Neo(anc)','Spec','Sub')


dnds_plot <- ggplot(all_dnds, aes(x=factor(category, levels = func_levels),y=MG94xREV_omega))+ 
  # facet_grid(cols=vars(variable)) +
  ylab('Dn / Ds') +
  geom_boxplot(color=c('#778899','#D4E056','#F67280',"#E9AB17","#E9AB17","#3B7BBD",'#E23F51','#E23F51','#6CBC4D','#F18244')) + 
  theme(axis.title.x=element_blank(),panel.background = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 20),
        strip.background =element_rect(color='black',fill="white")) +
  ylim(0.000000000000000000000000001,2) +
  #geom_jitter(size=0.001) + 
  stat_compare_means(aes(label = after_stat(p.signif)),method = "t.test", 
                     ref.group = 'All_Ortho',hide.ns = TRUE, size=5) +
  geom_hline(yintercept=median((all_dnds[all_dnds$category=='All_Ortho',])['MG94xREV_omega']$MG94xREV_omega),
             linetype="dashed", color = '#778899',linewidth=0.5) 


dn_plot <- ggplot(all_dnds, aes(x=factor(category, levels = func_levels),y=full_adaptive_model_Dn))+ 
  # facet_grid(cols=vars(variable)) +
  ylab('Dn') +
  geom_boxplot(color=c('#778899','#D4E056','#F67280',"#E9AB17","#E9AB17","#3B7BBD",'#E23F51','#E23F51','#6CBC4D','#F18244')) + 
  theme(axis.title.x=element_blank(),panel.background = element_rect(colour = "black", fill=NA, size=1),
        axis.text.x=element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 20),
        strip.background =element_rect(color='black',fill="white")) +
  ylim(0.000000000000000000000000001,2)  + 
  stat_compare_means(aes(label = after_stat(p.signif)),method = "t.test", 
                     ref.group = 'All_Ortho',hide.ns = TRUE, size=5) +
  geom_hline(yintercept=median((all_dnds[all_dnds$category=='All_Ortho',])['full_adaptive_model_Dn']$full_adaptive_model_Dn),
             linetype="dashed", color = '#778899',linewidth=0.5) 


ds_plot <- ggplot(all_dnds, aes(x=factor(category, levels = func_levels),y=full_adaptive_model_Ds))+ 
  # facet_grid(cols=vars(variable)) +
  ylab('Ds') +
  geom_boxplot(color=c('#778899','#D4E056','#F67280',"#E9AB17","#E9AB17","#3B7BBD",'#E23F51','#E23F51','#6CBC4D','#F18244')) + 
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

ggsave(filename = './Plots/dnds_dn_ds_boxplots.jpg', width = 7, height = 13)


# tau boxplot
all_tau <- merge(tau,all_dnds,by.x='YOgnID',by.y='id')


ggplot(all_tau, aes(x=factor(category, levels = func_levels),y=tau))+ 
  # facet_grid(cols=vars(variable)) +
  ylab('Tau') +
  geom_boxplot(color=c('#778899','#D4E056','#F67280',"#E9AB17","#E9AB17","#3B7BBD",'#E23F51','#E23F51','#6CBC4D','#F18244')) + 
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

ggsave(filename = './Plots/tau_boxplot.jpg', width = 7, height = 4)


# tissue expression heatmap 

tissue_exp <- pivot_longer(exp, cols = c('f_ac':'m_wb'))

all_tissue_exp <- merge(tissue_exp,all_tau,by='YOgnID')

all_tissue_exp <- all_tissue_exp %>% 
  select(YOgnID, name, value, category) %>% 
  group_by(name,category) %>%
  mutate(value=mean(value)) %>%
  select(name,category,value) %>%
  ungroup() %>%
  distinct()



ggplot(all_tissue_exp, aes(x=factor(category,levels=func_levels), y=name, fill= value)) +
  geom_tile() +
  scale_fill_gradient(high = "white", low = "dodgerblue4",guide = "colorbar") +
  theme(legend.title = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text=element_text(size=14),
        panel.background = element_rect(fill = "white", color = "white"))  

ggsave(filename = './Plots/expression_heatmap.jpg', width = 14, height = 7)



# p-value only analyses
p_val_analyses <- rbind(all_dups[c('id','category')],
                        all_ortho[c('id','category')])

# add exon amounts
colnames(exon_counts) <- c('id','n_exons')
p_val_analyses <- merge(p_val_analyses,exon_counts,by='id')

# add GC content and aa length 
seqs <- seqs[c('YOgn','prot_length','gc_content')]
colnames(seqs)[1] <- 'id'

p_val_analyses <- merge(p_val_analyses,seqs,by='id')

#


# absolute expression level
# sequence similarity ?





# correlation matrix plot 




##
# add GC content

# read in duplicate and ortho information
dups <- read.csv("C:/Users/17735/Downloads/Duplicate_Pairs_Project/Final_Figures/Input_Files/Dups_Master_File.tsv", sep="")
orthos <- read.csv("C:/Users/17735/Downloads/Duplicate_Pairs_Project/Final_Figures/Input_Files/Orthos_Master_File.tsv", sep="")

# read in fasta and exon information for all genes 
eight_species <- read.csv("C:/Users/17735/Downloads/Duplicate_Pairs_Project/Identifying_Duplicate_Pairs/Input_Files/EightSpecies.tsv", sep="")

library(Biostrings)

gc <- eight_species[,c(5,9)]
gc <- na.omit(gc)

colnames(gc) <- c('gn','nuc')

gc$nuc <- trimws(gc$nuc)

gc <- gc[!grepl("o", gc$nuc), ]

library(stringr)
gc$gc_content <- str_count(gc$nuc, "[GC]") / nchar(gc$nuc) * 100

gc <- gc[!duplicated(gc),]

# merge with dup and ortho dfs

colnames(gc)[1] <- 'dup_gn'
dups <- merge(dups,gc,by='dup_gn')

colnames(gc)[1] <- 'ortho_gn'
orthos <- merge(orthos,gc,by='ortho_gn')

dups_orig <- dups
orthos_orig <- orthos

## 
# add expression divergence
dups <- read.csv("C:/Users/17735/Downloads/Duplicate_Pairs_Project/Identifying_Duplicate_Pairs/Input_Files/Duplicate_Pairs_Information.tsv", sep="")
dups <- dups[,c(9,15)]

dup_1 <- dups[1]
dup_2 <- dups[2]

colnames(dup_1) <- 'gn'
colnames(dup_2) <- 'gn'

# read in expression data
exp <- read.csv("C:/Users/17735/Downloads/Duplicate_Pairs_Project/Identifying_Duplicate_Pairs/Input_Files/Expression_Data.tsv", sep="")

# merge duplicate copies with their expression data
colnames(exp)[1] <- 'gn'
dup_1 <- merge(dup_1,exp,by='gn')
dup_2 <- merge(dup_2,exp,by='gn')

# calculate euclidean distances between copies 
rownames(dup_1) <- dup_1$gn
dup_1 <- dup_1[,-c(1)]

rownames(dup_2) <- dup_2$gn
dup_2 <- dup_2[,-c(1)]

dup_ed <- (rowSums((dup_1 - dup_2) ^ 2)) ^ (1/2)
dups$ed <- dup_ed

# merge back 
#dups_orig <- read.csv("C:/Users/17735/Downloads/Duplicate_Pairs_Project/Final_Figures/Input_Files/Dups_Master_File.tsv", sep="")

dup1 <- dups[,c(1,3)]
dups <- dups[,c(2,3)]

colnames(dup1) <- c('dup_gn','Expression Divergence')
colnames(dups) <- c('dup_gn','Expression Divergence')

dups <- rbind(dups,dup1)

dups <- merge(dups_orig,dups,by='dup_gn')
##

# make correlation matrix for duplicates

dups_orig <- dups

library(ggcorrplot)
dups <- dups[,c(2:7,10:13,30,31)]
dups <- dups[!dups$dup_dnds==Inf,]


dup_cor <- cor(dups)
colnames(dup_cor) <- c('Exons','Length (bp)','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','GC Content','Expression Divergence')
rownames(dup_cor) <- c('Exons','Length (bp)','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','GC Content','Expression Divergence')
dup_cor <- dup_cor[, c('Exons','Length (bp)','GC Content','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','Expression Divergence')]
dup_cor <- dup_cor[c('Exons','Length (bp)','GC Content','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','Expression Divergence'),]
dup_pmat <- cor_pmat(dup_cor)

library(grDevices)
library(corrplot)
pdf(file="./Input_Files/Dup_Correlation_Matrix.pdf", height=10, width=10)
corrplot(dup_cor, p.mat=dup_pmat, method="circle", type="lower",sig.level=c(0.0001,0.001,0.01, 0.05),
         insig="label_sig", pch.col="white", tl.col="black",
         diag=FALSE, pch.cex=1.6, col=colorRampPalette(c("darkred","white","darkgreen"))(200))
dev.off()


##
dups <- dups_orig
dups <- dups[dups$func_copy == 'conserv', ]
dups <- dups[,c(2:7,10:13,30,31)]
dups <- dups[!dups$dup_dnds==Inf,]
dup_cor <- cor(dups)
colnames(dup_cor) <- c('Exons','Length (bp)','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','GC Content','Expression Divergence')
rownames(dup_cor) <- c('Exons','Length (bp)','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','GC Content','Expression Divergence')
dup_cor <- dup_cor[, c('Exons','Length (bp)','GC Content','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','Expression Divergence')]
dup_cor <- dup_cor[c('Exons','Length (bp)','GC Content','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','Expression Divergence'),]
dup_pmat <- cor_pmat(dup_cor)

pdf(file="./Input_Files/Conserved_Correlation_Matrix.pdf", height=10, width=10)
corrplot(dup_cor, p.mat=dup_pmat, method="circle", type="lower",sig.level=c(0.0001,0.001,0.01, 0.05),
         insig="label_sig", pch.col="white", tl.col="black",
         diag=FALSE, pch.cex=1.6, col=colorRampPalette(c("darkred","white","darkgreen"))(200))
dev.off()

##
dups <- dups_orig
dups <- dups[dups$func_copy == 'neo_copy',]
dups <- dups[,c(2:7,10:13,30,31)]
dups <- dups[!dups$dup_dnds==Inf,]
dup_cor <- cor(dups)
colnames(dup_cor) <- c('Exons','Length (bp)','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','GC Content','Expression Divergence')
rownames(dup_cor) <- c('Exons','Length (bp)','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','GC Content','Expression Divergence')
dup_cor <- dup_cor[, c('Exons','Length (bp)','GC Content','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','Expression Divergence')]
dup_cor <- dup_cor[c('Exons','Length (bp)','GC Content','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','Expression Divergence'),]
dup_pmat <- cor_pmat(dup_cor)

pdf(file="./Input_Files/NeoCopy_Correlation_Matrix.pdf", height=10, width=10)
corrplot(dup_cor, p.mat=dup_pmat, method="circle", type="lower",sig.level=c(0.0001,0.001,0.01, 0.05),
         insig="label_sig", pch.col="white", tl.col="black",
         diag=FALSE, pch.cex=1.6, col=colorRampPalette(c("darkred","white","darkgreen"))(200))
dev.off()

##
dups <- dups_orig
dups <- dups[dups$func_copy == 'not_neo_copy',]
dups <- dups[,c(2:7,10:13,30,31)]
dups <- dups[!dups$dup_dnds==Inf,]
dup_cor <- cor(dups)
colnames(dup_cor) <- c('Exons','Length (bp)','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','GC Content','Expression Divergence')
rownames(dup_cor) <- c('Exons','Length (bp)','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','GC Content','Expression Divergence')
dup_cor <- dup_cor[, c('Exons','Length (bp)','GC Content','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','Expression Divergence')]
dup_cor <- dup_cor[c('Exons','Length (bp)','GC Content','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','Expression Divergence'),]
dup_pmat <- cor_pmat(dup_cor)

pdf(file="./Input_Files/NotNeoCopy_Correlation_Matrix.pdf", height=10, width=10)
corrplot(dup_cor, p.mat=dup_pmat, method="circle", type="lower",sig.level=c(0.0001,0.001,0.01, 0.05),
         insig="label_sig", pch.col="white", tl.col="black",
         diag=FALSE, pch.cex=1.6, col=colorRampPalette(c("darkred","white","darkgreen"))(200))
dev.off()

##
dups <- dups_orig
dups <- dups[dups$func_copy == 'specializ',]
dups <- dups[,c(2:7,10:13,30,31)]
dups <- dups[!dups$dup_dnds==Inf,]
dup_cor <- cor(dups)
colnames(dup_cor) <- c('Exons','Length (bp)','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','GC Content','Expression Divergence')
rownames(dup_cor) <- c('Exons','Length (bp)','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','GC Content','Expression Divergence')
dup_cor <- dup_cor[, c('Exons','Length (bp)','GC Content','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','Expression Divergence')]
dup_cor <- dup_cor[c('Exons','Length (bp)','GC Content','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','Expression Divergence'),]
dup_pmat <- cor_pmat(dup_cor)

pdf(file="./Input_Files/Specializ_Correlation_Matrix.pdf", height=10, width=10)
corrplot(dup_cor, p.mat=dup_pmat, method="circle", type="lower",sig.level=c(0.0001,0.001,0.01, 0.05),
         insig="label_sig", pch.col="white", tl.col="black",
         diag=FALSE, pch.cex=1.6, col=colorRampPalette(c("darkred","white","darkgreen"))(200))
dev.off()




#####################################################################################


## repeat for orthologs
#orthos <- orthos[,c(2:16)]

# add expression divergence
orthos <- read.csv("C:/Users/17735/Downloads/Duplicate_Pairs_Project/Final_Figures/Input_Files/Orthos_Master_File.tsv", sep="")


#orthos <- orthos[,c(9,15)]

ortho_1 <- orthos[,c(1,17,18)]
ortho_2 <- ortho_1[ortho_1$ortho=='ortho_2',]
ortho_1 <- ortho_1[ortho_1$ortho=='ortho_1',]

ortho_1 <- ortho_1[1]
ortho_2 <- ortho_2[1]

colnames(ortho_1) <- 'gn'
colnames(ortho_2) <- 'gn'

# merge ortholicate copies with their expression data
ortho_1 <- merge(ortho_1,exp,by='gn')
ortho_2 <- merge(ortho_2,exp,by='gn')

# calculate euclidean distances between copies 
rownames(ortho_1) <- ortho_1$gn
ortho_1 <- ortho_1[,-c(1)]

rownames(ortho_2) <- ortho_2$gn
ortho_2 <- ortho_2[,-c(1)]

ortho_ed <- (rowSums((ortho_1 - ortho_2) ^ 2)) ^ (1/2)
orthos$ed <- ortho_ed

# merge back 
#orthos_orig <- read.csv("C:/Users/17735/Downloads/Duplicate_Pairs_Project/Final_Figures/Input_Files/Orthos_Master_File.tsv", sep="")

ortho1 <- orthos[,c(1,3)]
orthos <- orthos[,c(2,3)]

colnames(ortho1) <- c('ortho_gn','Expression Divergence')
colnames(orthos) <- c('ortho_gn','Expression Divergence')

orthos <- rbind(orthos,ortho1)

orthos <- merge(orthos_orig,orthos,by='ortho_gn')
##

# make correlation matrix for orthologs

library(ggcorrplot)
orthos <- orthos[,c(2:11,20,21)]
orthos <- orthos[!orthos$ortho_dnds==Inf,]



ortho_cor <- cor(orthos)

colnames(ortho_cor) <- c('Tau','Tissues Expressed in','Expression Level','Exons','Length (bp)','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','GC Content','Expression Divergence')
rownames(ortho_cor) <- c('Tau','Tissues Expressed in','Expression Level','Exons','Length (bp)','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','GC Content','Expression Divergence')
ortho_cor <- ortho_cor[, c('Exons','Length (bp)','GC Content','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','Expression Divergence')]
ortho_cor <- ortho_cor[c('Exons','Length (bp)','GC Content','Sequence Similarity','Dn','Ds','Dn/Ds','PPI','Tau','Tissues Expressed in','Expression Level','Expression Divergence'),]

ortho_pmat <- cor_pmat(ortho_cor)


library(grDevices)
library(corrplot)
pdf(file="./Input_Files/Ortho_Correlation_Matrix.pdf", height=10, width=10)
corrplot(ortho_cor, p.mat=ortho_pmat, method="circle", type="lower",sig.level=c(0.0001,0.001,0.01, 0.05),
         insig="label_sig", pch.col="white", tl.col="black",
         diag=FALSE, pch.cex=1.6, col=colorRampPalette(c("darkred","white","darkgreen"))(200))
dev.off()




# n orthologs before pval analyses
# N ORTHOLOGS AND DS AND FUNC 













