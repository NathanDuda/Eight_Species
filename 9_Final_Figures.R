


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

# calculate p values
n_exons_pval <- t.test(p_val_analyses$n_exons[p_val_analyses$category=='All_Dup'],
                       p_val_analyses$n_exons[p_val_analyses$category=='All_Ortho'])
n_exons_pval <- n_exons_pval$p.value
n_exons_dup <- mean(p_val_analyses$n_exons[p_val_analyses$category=='All_Dup'])
n_exons_ortho <- mean(p_val_analyses$n_exons[p_val_analyses$category=='All_Ortho'])


gc_content_pval <- t.test(p_val_analyses$gc_content[p_val_analyses$category=='All_Dup'],
                          p_val_analyses$gc_content[p_val_analyses$category=='All_Ortho'])
gc_content_pval <- gc_content_pval$p.value
gc_content_dup <- mean(p_val_analyses$gc_content[p_val_analyses$category=='All_Dup'])
gc_content_ortho <- mean(p_val_analyses$gc_content[p_val_analyses$category=='All_Ortho'])


prot_length_pval <- t.test(p_val_analyses$prot_length[p_val_analyses$category=='All_Dup'],
                           p_val_analyses$prot_length[p_val_analyses$category=='All_Ortho'])
prot_length_pval <- prot_length_pval$p.value
prot_length_dup <- mean(p_val_analyses$prot_length[p_val_analyses$category=='All_Dup'])
prot_length_ortho <- mean(p_val_analyses$prot_length[p_val_analyses$category=='All_Ortho'])


# absolute mean expression level 
tissue_exp <- exp %>%
  pivot_longer(., cols = c('f_ac':'m_wb')) %>%
  group_by(YOgnID) %>%
  mutate(value = sum(value)) %>%
  select(-name) %>%
  distinct()

colnames(tissue_exp)[1] <- 'id'
abs_exp <- merge(p_val_analyses, tissue_exp, by = 'id')


abs_exp_pval <- t.test(abs_exp$value[abs_exp$category=='All_Dup'],
                       abs_exp$value[abs_exp$category=='All_Ortho'])
abs_exp_pval <- abs_exp_pval$p.value
abs_exp_dup <- mean(abs_exp$value[abs_exp$category=='All_Dup'])
abs_exp_ortho <- mean(abs_exp$value[abs_exp$category=='All_Ortho'])


### absolute expression ratio of duplicate pairs versus ancestral copy

# get expression values for each duplicate copy
dup1_exp <- exp %>% rename_at(-1, ~paste('dup1_', ., sep = '')) %>% rename_with(~ paste0("dup_1", names(.)[1]), 1)
dups_exp <- merge(dup1_exp,dups,by='dup_1')

dup2_exp <- exp %>% rename_at(-1, ~paste('dup2_', ., sep = '')) %>% rename_with(~ paste0("dup_2", names(.)[1]), 1)
dups_exp <- merge(dup2_exp,dups_exp,by='dup_2')

anc_exp <- exp %>% rename_at(-1, ~paste('anc_', ., sep = '')) %>% rename_with(~ paste0("ancestral_copy", names(.)[1]), 1)
dups_exp <- merge(anc_exp,dups_exp,by='ancestral_copy')


# calculate absolute expression ratio 
ratio <- dups_exp %>%
  select(-dup_1, -dup_2, -ancestral_copy) %>%
  mutate(dup1_sum = rowSums(select(., starts_with('dup1')))) %>%
  mutate(dup2_sum = rowSums(select(., starts_with('dup2')))) %>%
  mutate(anc_sum = rowSums(select(., starts_with('anc_')))) %>%
  select(Orthogroup, dup1_sum, dup2_sum, anc_sum) %>%
  mutate(ratio=(dup1_sum + dup2_sum) / anc_sum) 

median(ratio$ratio)
# 1:1.05 ratio compared to 1:2, suggesting dosage is in effect (explains conservation)

sd(ratio$ratio)
# but huge standard deviation 

ggplot(ratio, aes(x=ratio)) +
  geom_density() +
  xlim(0,10) +
  geom_vline(xintercept =1)

ratio <- ratio[c('Orthogroup','ratio')]


# ratio per tissue
tissue_ratio <- dups_exp %>%
  select(-dup_1, -dup_2, -ancestral_copy) 

tissue_names <- c('f_ac','f_dg','f_go','f_hd','f_re','f_tx','f_wb',
                  'm_ac','m_dg','m_go','m_hd','m_re','m_tx','m_wb')

for (tissue in tissue_names) {
  tissue_ratio <- tissue_ratio %>%
    mutate(!!tissue := (select(., paste0('dup1_', !!tissue))   +
                          select(., paste0('dup2_', !!tissue))  ) / 
                          select(., paste0('anc_', !!tissue)) ) %>%
    select(-paste0('dup1_', tissue),
           -paste0('dup2_', tissue),
           -paste0('anc_', tissue)) %>%
    mutate_all(~unlist(.))
}

tissue_ratio <- tissue_ratio %>%
  select(-Orthogroup) %>%
  mutate_all(~ replace(., is.na(.), NA) %>% 
               replace(., is.infinite(.), NA)) %>%
  summarise_all(~ median(., na.rm = TRUE)) %>%
  t()

###


# get number of genes in each orthogroup
n_gene_per_orthogroup <- orthogroups %>%
  mutate_all(~na_if(., "")) %>%
  mutate(n_genes = 9 - rowSums(is.na(.))) %>%
  select(Orthogroup,n_genes)

# correlation matrix plot 
corr_plot <- merge(tissue_exp,p_val_analyses, by = 'id')
corr_plot <- merge(corr_plot,tau,by.x='id', by.y='YOgnID')

dup_corr_plot <- merge(corr_plot,all_dups,by='id')
dup_corr_plot <- merge(dup_corr_plot,n_gene_per_orthogroup,by='Orthogroup')
dup_corr_plot <- merge(dup_corr_plot,ratio, by = 'Orthogroup')

ortho_corr_plot <- merge(corr_plot,all_ortho,by='id')
ortho_corr_plot <- merge(ortho_corr_plot,n_gene_per_orthogroup,by='Orthogroup')
ortho_corr_plot <- merge(ortho_corr_plot,ratio, by = 'Orthogroup')


# make correlation matrix for duplicates
dup_corr_plot <- dup_corr_plot[c('n_exons','prot_length','full_adaptive_model_Dn','full_adaptive_model_Ds','MG94xREV_omega','tau','value','gc_content','n_genes','ratio')]
dup_corr_plot <- dup_corr_plot %>% mutate_all(~ as.numeric(.))

library(ggcorrplot)
dup_corr_plot <- cor(dup_corr_plot)
colnames(dup_corr_plot) <- c('Exon Amount','Length (bp)','Dn','Ds','Dn/Ds','Tau','Expression Level','GC Content','Genes in Orthogroup','Expression Ratio')
rownames(dup_corr_plot) <- colnames(dup_corr_plot)

dup_pmat <- cor_pmat(dup_corr_plot)

library(grDevices)
library(corrplot)
pdf(file="./Plots/dup_correlation_matrix.pdf", height=5, width=5)
corrplot(dup_corr_plot, p.mat=dup_pmat, method="circle", type="lower",sig.level=c(0.0001,0.001,0.01, 0.05),
         insig="label_sig", pch.col="white", tl.col="black",
         diag=FALSE, pch.cex=1.6, col=colorRampPalette(c("darkred","white","darkgreen"))(200))
dev.off()


# make correlation matrix for orthologs
ortho_corr_plot <- ortho_corr_plot[c('n_exons','prot_length','full_adaptive_model_Dn','full_adaptive_model_Ds','MG94xREV_omega','tau','value','gc_content','n_genes','ratio')]
ortho_corr_plot <- ortho_corr_plot %>% mutate_all(~ as.numeric(.))

ortho_corr_plot <- cor(ortho_corr_plot)
colnames(ortho_corr_plot) <- c('Exon Amount','Length (bp)','Dn','Ds','Dn/Ds','Tau','Expression Level','GC Content','Genes in Orthogroup','Expression Ratio')
rownames(ortho_corr_plot) <- colnames(ortho_corr_plot)

ortho_pmat <- cor_pmat(ortho_corr_plot)

pdf(file="./Plots/ortho_correlation_matrix.pdf", height=5, width=5)
corrplot(ortho_corr_plot, p.mat=ortho_pmat, method="circle", type="lower",sig.level=c(0.0001,0.001,0.01, 0.05),
         insig="label_sig", pch.col="white", tl.col="black",
         diag=FALSE, pch.cex=1.6, col=colorRampPalette(c("darkred","white","darkgreen"))(200))
dev.off()



# phylogeny with ds histogram
ds_tree <- dups_dnds[c('id','full_adaptive_model_Ds')]
ds_tree <- merge(ds_tree,longest_transcript[c('YOgn','species')],by.x='id',by.y='YOgn')

tree <- ape::read.tree('./MEGA_Tree.nwk')

library(ggtree)
ggplot(ds_tree, aes(x=full_adaptive_model_Ds,y=factor(species, levels=rev(c('dmel','dyak','dana','dpse','dper','dvir','dmoj','dwil'))))) +
  geom_jitter() +
  xlim(2,0) +
  ylab('') +
  xlab('Ds') +
  theme_bw() +
  theme(panel.border = element_blank()) +
  ggtree(tree) + geom_tiplab(aes(label = label)) 

ggsave(filename = './Plots/ds_phylogeny.jpg', width = 13, height = 5)




