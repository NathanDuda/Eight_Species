


source('./startup.R')


dups <- read.csv("./Dup_Pairs_Ancestral.tsv", sep = '')

orthogroups <- read.csv("./Dup_Pair_Orthologs.tsv", sep = '')

exp_rel <- read.csv("./Relative_Expression.tsv", sep = '')
exp <- read.csv("./All_Expression_Data.tsv", sep = '')

tau <- read.csv("./Tau.tsv", sep = '')

dnds <- read.csv("./orthogroups_omega.tsv", sep = '')

longest_transcript <- read.csv("./longest_transcript.tsv", sep = '')
seqs <- read.csv("./longest_transcript.tsv", sep = '')
seqs <- seqs %>% mutate(gc_content = (str_count(longest_ORF, "G") + str_count(longest_ORF, "C")) / nchar(longest_ORF) * 100)

exon_counts <- read.csv("./Exon_counts.tsv", sep = '')

mech <- read.csv("./Dup_Mechanism.tsv", sep = "")
mech <- mech[c('Orthogroup','mech')]






dups_orthogroups <- read.csv("./Duplicate_Pairs.tsv", sep="")


func <- read.csv("./CDROM_CLOUD_pseudo_funcs.tsv", sep = "") %>%
  select(dup_1, dup_2, anc, func_pseudo = CLOUD_pseudo) %>% # Change to CLOUD when desired
  merge(., dups_orthogroups, by = c('dup_1', 'dup_2'))


orig_func <- func
############################


func <- func %>% select(Orthogroup, func_pseudo)
func$func_pseudo <- gsub('conserv', 'conserved', func$func_pseudo)
func$func_pseudo <- gsub('specializ', 'specialized', func$func_pseudo)
func$func_pseudo <- gsub('subfun', 'subfunctionalized', func$func_pseudo)
func$func_pseudo <- gsub('neo', 'neofunctionalized', func$func_pseudo)


# per orthogroup figures:


# func x mech table 
orthogroup <- merge(dups, func, by = 'Orthogroup')
orthogroup <- merge(orthogroup, mech, by = 'Orthogroup')

orthogroup <- orthogroup[orthogroup$Orthogroup %in% dnds$Orthogroup,]

func_mech_table <- as.data.frame(unclass(table(orthogroup$mech, orthogroup$func)))
colnames(func_mech_table) <- c('conserved','neo_dup1','neo_dup2','specialized','subfunctionalized')

library(knitr)
library(kableExtra)
func_mech_table %>% 
  kable(align = rep('c', 5)) %>%
  column_spec (1:6,border_left = T, border_right = T) %>%
  kable_styling()


# func pie 
func_pie <- orthogroup %>% select(func = func_pseudo)

func_pie$func <- gsub('neofunctionalized_dup1', 'neofunctionalized', func_pie$func)
func_pie$func <- gsub('neofunctionalized_dup2', 'neofunctionalized', func_pie$func)
func_pie$func <- gsub('pseudo_dup_1', 'pseudo', func_pie$func)
func_pie$func <- gsub('pseudo_dup_2', 'pseudo', func_pie$func)

func_pie <- as.data.frame(table(func_pie$func))

colnames(func_pie) <- c('func', 'value')

ggplot(func_pie, aes(x = "", y = value, group = func, fill = func)) +
  geom_bar(width = 1, stat = "identity", position = position_fill()) +
  geom_text(aes(label = value), position = position_fill(vjust = 0.5), colour = 'black', size = 3) +
  theme_void() +
  theme(legend.title = element_blank(), panel.background = element_rect(fill = "white", color = "white"), legend.background = element_rect(fill = "white", color = 'white')) +
  coord_polar("y") +
  scale_fill_manual(values=c("#3B7BBD", "#E23F51", 'gray', 'darkgray', "#6CBC4D", '#F18244'))

ggsave(filename = './Plots/func_pie.jpg', width = 5, height = 3, device='tiff', dpi=300)


# mech pie 
mech_pie <- orthogroup %>% select(mech)

mech_pie$mech <- gsub('_dup_1', '', mech_pie$mech)
mech_pie$mech <- gsub('_dup_2', '', mech_pie$mech)

mech_pie <- as.data.frame(table(mech_pie$mech))
colnames(mech_pie) <- c('mech','value')

ggplot(mech_pie, aes(x = "", y = value, group = mech, fill = mech)) +
  geom_bar(width = 1, stat = "identity", position = position_fill()) +
  geom_text(aes(label = value), position = position_fill(vjust = 0.5), colour = 'black', size = 3) +
  theme_void() +
  theme(legend.title = element_blank(), panel.background = element_rect(fill = "white", color = "white"), legend.background = element_rect(fill = "white", color = 'white')) +
  coord_polar("y") +
  scale_fill_manual(values=c("#F67280", "#E9AB17", "#1E90FF"))

ggsave(filename = './Plots/mech_pie.jpg', width = 5, height = 3, device='tiff', dpi=300)


# ds group func barplot 
dups_ds_func <- dnds %>%
  mutate(Ds_group = case_when(FAM_ds <= 0.125 & FAM_ds > 0 ~ '0.125',
                              FAM_ds <= 0.25 & FAM_ds > 0.125 ~ '0.25',
                              FAM_ds <= 0.375 & FAM_ds > 0.25 ~ '0.375',
                              FAM_ds <= 0.5 & FAM_ds > 0.375 ~ '0.5',
                              FAM_ds <= 0.625 & FAM_ds > 0.5 ~ '0.625',
                              FAM_ds <= 0.75 & FAM_ds > 0.625 ~ '0.75',
                              FAM_ds <= 0.875 & FAM_ds > 0.75 ~ '0.875',
                              FAM_ds <= 1 & FAM_ds > 0.875 ~ '1',
                              FAM_ds > 1 ~ '>1',))

dups_ds_func <- merge(func, dups_ds_func, by='Orthogroup')

dups_ds_func$func <- gsub('neofunctionalized_dup1', 'neofunctionalized', dups_ds_func$func)
dups_ds_func$func <- gsub('neofunctionalized_dup2', 'neofunctionalized', dups_ds_func$func)
dups_ds_func$func <- gsub('pseudo_dup_1', 'pseudo', dups_ds_func$func)
dups_ds_func$func <- gsub('pseudo_dup_2', 'pseudo', dups_ds_func$func)

dups_ds_func <- as.data.frame(table(dups_ds_func$func, dups_ds_func$Ds_group))

# filled barplot
ggplot(dups_ds_func, aes(fill = Var1, y = Freq, x = factor(Var2, levels = c('0.125','0.25','0.375','0.5','0.625','0.75','0.875','1','>1')))) + 
  geom_bar(position = "fill", stat = "identity") +
  scale_fill_manual(values = c("#3B7BBD", '#E23F51', 'gray', 'darkgray', '#6CBC4D', '#F18244'))+
  ylab('Percentage') +
  xlab('Ds') +
  theme_classic() +
  theme(legend.title = element_blank())

ggsave(filename = './Plots/ds_grouped_func_filled_barplot.jpg', width = 7, height = 4, device='tiff', dpi=300)

# stacked barplot
ggplot(dups_ds_func, aes(fill = Var1, y = Freq, x = factor(Var2, levels = c('0.125','0.25','0.375','0.5','0.625','0.75','0.875','1','>1')))) + 
  geom_bar(position = "stack", stat = "identity") +
  scale_fill_manual(values = c("#3B7BBD", '#E23F51', 'gray', 'darkgray', '#6CBC4D', '#F18244'))+
  ylab('Number') +
  xlab('Ds') +
  theme_classic() +
  theme(legend.title = element_blank())

ggsave(filename = './Plots/ds_grouped_func_stacked_barplot.jpg', width = 7, height = 4, device='tiff', dpi=300)


######################################################################################
# dup copy figures: 
dups_dnds <- dnds %>% 
  select(Orthogroup, id, baseline_MG94xREV_omega, FAM_dn, FAM_ds, dup_or_ortho)
  
# combine dn, ds, dnds, func, and mech into one dataframe

dups_dnds <- full_join(dups_dnds, func, by = 'Orthogroup')
dups_dnds <- full_join(dups_dnds, mech, by = 'Orthogroup')
dups_dnds <- full_join(dups_dnds, dups, by = 'Orthogroup')
# some one to two orthogroups not matching to dup pair bc ancestral copy not expressed 
# some dup pairs not matching to orthogroups because they dont have at least 4 genes (required by hyphy)


# separate duplicate genes and ortholog genes in the orthogroups into their own two dataframes 
all_dups <- dups_dnds[dups_dnds$dup_or_ortho == 'dup' & !is.na(dups_dnds$dup_or_ortho),]
all_dups$category <- 'All_Dup'
all_ortho <- dups_dnds[dups_dnds$dup_or_ortho == 'ortho' & !is.na(dups_dnds$dup_or_ortho),]
all_ortho$category <- 'All_Ortho'

all_func <- dups_dnds %>%
  filter(!is.na(func_pseudo) & dup_or_ortho == 'dup') %>%
  mutate(category = case_when(func_pseudo == 'neofunctionalized_dup1' & dup_1 == id ~ 'Neo(new)',
                              func_pseudo == 'neofunctionalized_dup1' & dup_2 == id ~ 'Neo(anc)',
                              func_pseudo == 'neofunctionalized_dup2' & dup_1 == id ~ 'Neo(anc)',
                              func_pseudo == 'neofunctionalized_dup2' & dup_2 == id ~ 'Neo(new)',
                              func_pseudo == 'cons' ~ 'Cons',
                              func_pseudo == 'spec' ~ 'Spec',
                              func_pseudo == 'sub' ~ 'Sub',
                              func_pseudo == 'pseudo_dup_1' & dup_1 == id ~ 'Pseudo(non)',
                              func_pseudo == 'pseudo_dup_1' & dup_2 == id ~ 'Pseudo(anc)',
                              func_pseudo == 'pseudo_dup_2' & dup_1 == id ~ 'Pseudo(anc)',
                              func_pseudo == 'pseudo_dup_2' & dup_2 == id ~ 'Pseudo(non)',
                              func_pseudo == 'pseudo_both' ~ 'Pseudo(both)'))

all_mech <- dups_dnds %>%
  filter(!is.na(mech) & dup_or_ortho =='dup' & mech !='unknown') %>%
  mutate(category = case_when(mech == 'RNA_dup_1' & dup_1 == id ~ 'RNA(new)',
                              mech == 'RNA_dup_1' & dup_2 == id ~ 'RNA(anc)',
                              mech == 'RNA_dup_2' & dup_1 == id ~ 'RNA(anc)',
                              mech == 'RNA_dup_2' & dup_2 == id ~ 'RNA(new)',
                              T ~ mech))

# make one master dataframe with all the categories
all_dnds <- rbind(all_dups, all_ortho, all_func, all_mech)


func_levels <- c('All_Ortho','All_Dup','DNA','RNA(new)','RNA(anc)','Cons','Neo(new)',
                 'Neo(anc)','Pseudo(non)','Pseudo(anc)','Pseudo(both)','Spec','Sub')


dnds_plot <- all_dnds %>%
  filter(category %in% func_levels) %>%
  ggplot(., aes(x = factor(category, levels = func_levels), y = baseline_MG94xREV_omega))+ 
  # facet_grid(cols=vars(variable)) +
  ylab('Dn / Ds') +
    geom_boxplot(color=c('black','#D4E056','#0079C2',"#FF828B","#FF828B","#778899",'#C05780','#C05780','lightgray','gray','darkgray','#E9AB17','#6CBC4D')) + 
  theme(axis.title.x=element_blank(),panel.background = element_rect(colour = "black", fill=NA, size=1),
       # axis.text.x=element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 20),
        strip.background =element_rect(color='black',fill="white")) +
  ylim(0.000000000000000000000000001,2) +
  #geom_jitter(size=0.001) + 
  stat_compare_means(aes(label = after_stat(p.signif)),method = "t.test", 
                     ref.group = 'All_Ortho',hide.ns = TRUE, size=5) +
  geom_hline(yintercept=median((all_dnds[all_dnds$category=='All_Ortho',])['baseline_MG94xREV_omega']$baseline_MG94xREV_omega),
             linetype="dashed", color = '#778899',linewidth=0.5) 

ggsave(dnds_plot, filename = './Plots/dnds_boxplot.jpg', width = 7, height = 4)


dn_plot <- ggplot(all_dnds, aes(x = factor(category, levels = func_levels), y = FAM_dn))+ 
  # facet_grid(cols=vars(variable)) +
  ylab('Dn') +
  geom_boxplot(color = c('black','#D4E056','#0079C2',"#FF828B","#FF828B","#778899",'#C05780','#C05780','lightgray','gray','darkgray','#E9AB17','#6CBC4D')) + 
  theme(axis.title.x = element_blank(), panel.background = element_rect(colour = "black", fill = NA, size = 1),
       # axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 20),
        strip.background = element_rect(color = 'black', fill = "white")) +
  ylim(0.000000000000000000000000001, 2)  + 
  stat_compare_means(aes(label = after_stat(p.signif)), method = "t.test", 
                     ref.group = 'All_Ortho',hide.ns = TRUE, size = 5) +
  geom_hline(yintercept = median((all_dnds[all_dnds$category == 'All_Ortho',])['FAM_dn']$FAM_dn),
             linetype = "dashed", color = '#778899', linewidth = 0.5) 


ds_plot <- ggplot(all_dnds, aes(x = factor(category, levels = func_levels), y = FAM_ds))+ 
  # facet_grid(cols=vars(variable)) +
  ylab('Ds') +
  geom_boxplot(color = c('black','#D4E056','#0079C2',"#FF828B","#FF828B","#778899",'#C05780','#C05780','lightgray','gray','darkgray','#E9AB17','#6CBC4D')) + 
  theme(axis.title.x = element_blank(), panel.background = element_rect(colour = "black", fill = NA, size = 1),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 20),
        strip.background = element_rect(color = 'black', fill = "white")) +
  ylim(0.000000000000000000000000001, 2)  + 
  stat_compare_means(aes(label = after_stat(p.signif)), method = "t.test", 
                     ref.group = 'All_Ortho',hide.ns = TRUE, size = 5) +
  geom_hline(yintercept = median((all_dnds[all_dnds$category == 'All_Ortho',])['FAM_ds']$FAM_ds),
             linetype = "dashed", color = '#778899', linewidth = 0.5) 



ggarrange(dnds_plot, dn_plot, ds_plot, nrow = 3, ncol = 1)

ggsave(filename = './Plots/dnds_dn_ds_boxplots.jpg', width = 7, height = 13)


# tau boxplot
all_tau <- merge(tau, all_dnds, by.x = 'YOgnID', by.y = 'id')


ggplot(all_tau, aes(x = factor(category, levels = func_levels), y = tau))+ 
  # facet_grid(cols = vars(variable)) +
  ylab('Tau') +
  geom_boxplot(color = c('black','#D4E056','#0079C2',"#FF828B","#FF828B","#778899",'#C05780','#C05780','gray','#E9AB17','#6CBC4D')) + 
  theme(axis.title.x = element_blank(), panel.background = element_rect(colour = "black", fill = NA, size = 1),
        #axis.text.x = element_blank(),
        axis.text.y = element_text(size = 13),
        axis.title.y = element_text(size = 20),
        strip.background = element_rect(color = 'black', fill = "white")) +
  ylim(0, 1) +
  #geom_jitter(size = 0.001) + 
  stat_compare_means(aes(label = after_stat(p.signif)),method = "t.test", 
                     ref.group = 'All_Ortho',hide.ns = TRUE, size = 5) +
  geom_hline(yintercept = median((all_tau[all_tau$category == 'All_Ortho',])['tau']$tau),
             linetype = "dashed", color = '#778899',linewidth = 0.5) 
# pseudo(non) and pseudo(both) not included because have no expression 

ggsave(filename = './Plots/tau_boxplot.jpg', width = 7, height = 4)

write.table(all_tau, 'all_plotting_categories.tsv') # REMOVE IF PPI REMOVED OR MOVED INTO THIS SCRIPT


# tissue expression heatmap 

tissue_exp <- pivot_longer(exp_rel, cols = c('f_ac':'m_wb'))

all_tissue_exp <- merge(tissue_exp, all_tau, by = 'YOgnID')

all_tissue_exp <- all_tissue_exp %>% 
  select(YOgnID, name, value, category) %>% 
  group_by(name, category) %>%
  mutate(value = mean(value)) %>%
  select(name, category, value) %>%
  ungroup() %>%
  distinct()



ggplot(all_tissue_exp, aes(x = factor(category, levels = func_levels), y = name, fill = value)) +
  geom_tile() +
  scale_fill_gradient(high = "dodgerblue4", low = "#9AC0FF", guide = "colorbar") +
  theme(legend.title = element_blank(),
        strip.text = element_blank(),
        axis.title.x = element_blank(), axis.title.y = element_blank(),
        axis.text = element_text(size = 14),
        panel.background = element_rect(fill = "white", color = "white"))  

ggsave(filename = './Plots/expression_heatmap.jpg', width = 14, height = 7)



# p-value only analyses
p_val_analyses <- rbind(all_dups[c('id','category')],
                        all_ortho[c('id','category')])

# add exon amounts
colnames(exon_counts) <- c('id','n_exons')
p_val_analyses <- merge(p_val_analyses, exon_counts, by = 'id')

# add GC content and aa length 
seqs <- seqs %>% select(id = YOgn, prot_length, gc_content)

p_val_analyses <- merge(p_val_analyses, seqs, by = 'id')

# calculate p values
n_exons_pval <- t.test(p_val_analyses$n_exons[p_val_analyses$category == 'All_Dup'],
                       p_val_analyses$n_exons[p_val_analyses$category == 'All_Ortho'])
n_exons_pval <- n_exons_pval$p.value
n_exons_dup <- mean(p_val_analyses$n_exons[p_val_analyses$category == 'All_Dup'])
n_exons_ortho <- mean(p_val_analyses$n_exons[p_val_analyses$category == 'All_Ortho'])


gc_content_pval <- t.test(p_val_analyses$gc_content[p_val_analyses$category == 'All_Dup'],
                          p_val_analyses$gc_content[p_val_analyses$category == 'All_Ortho'])
gc_content_pval <- gc_content_pval$p.value
gc_content_dup <- mean(p_val_analyses$gc_content[p_val_analyses$category == 'All_Dup'])
gc_content_ortho <- mean(p_val_analyses$gc_content[p_val_analyses$category == 'All_Ortho'])


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
dups_exp <- merge(dup1_exp, dups, by = 'dup_1')

dup2_exp <- exp %>% rename_at(-1, ~paste('dup2_', ., sep = '')) %>% rename_with(~ paste0("dup_2", names(.)[1]), 1)
dups_exp <- merge(dup2_exp, dups_exp, by = 'dup_2')

anc_exp <- exp %>% rename_at(-1, ~paste('anc_', ., sep = '')) %>% rename_with(~ paste0("ancestral_copy", names(.)[1]), 1)
dups_exp <- merge(anc_exp, dups_exp, by = 'ancestral_copy')


# calculate absolute expression ratio 
ratio <- dups_exp %>%
  select(-dup_1, -dup_2, -ancestral_copy) %>%
  mutate(dup1_sum = rowSums(select(., starts_with('dup1')))) %>%
  mutate(dup2_sum = rowSums(select(., starts_with('dup2')))) %>%
  mutate(anc_sum = rowSums(select(., starts_with('anc_')))) %>%
  select(Orthogroup, dup1_sum, dup2_sum, anc_sum) %>%
  mutate(ratio = (dup1_sum + dup2_sum) / anc_sum) 

median(ratio$ratio)
# 1:0.9513249 ratio compared to 1:2, suggesting dosage is in effect (explains conservation)

sd(ratio$ratio)
# but huge standard deviation 

ggplot(ratio, aes(x = ratio)) +
  geom_density() +
  xlim(0, 10) +
  geom_vline(xintercept = 1)

ratio <- ratio[c('Orthogroup','ratio')]


# ratio per tissue
tissue_ratio <- dups_exp %>%
  select(-dup_1, -dup_2, -ancestral_copy, -pseudo) 

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
  select(Orthogroup, n_genes)

# correlation matrix plot 
corr_plot <- merge(tissue_exp, p_val_analyses, by = 'id')
corr_plot <- merge(corr_plot, tau, by.x = 'id', by.y = 'YOgnID')

dup_corr_plot <- merge(corr_plot, all_dups, by = 'id')
dup_corr_plot <- merge(dup_corr_plot, n_gene_per_orthogroup,by = 'Orthogroup')
dup_corr_plot <- merge(dup_corr_plot, ratio, by = 'Orthogroup')

ortho_corr_plot <- merge(corr_plot, all_ortho, by = 'id')
ortho_corr_plot <- merge(ortho_corr_plot, n_gene_per_orthogroup, by = 'Orthogroup')
ortho_corr_plot <- merge(ortho_corr_plot, ratio, by = 'Orthogroup')


# make correlation matrix for duplicates
dup_corr_plot <- dup_corr_plot %>% 
  select(n_exons, prot_length, FAM_dn, FAM_ds, baseline_MG94xREV_omega, tau, value, gc_content, n_genes, ratio) %>% 
  mutate_all(~ as.numeric(.))

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
ortho_corr_plot <- ortho_corr_plot %>%
  select(n_exons, prot_length, FAM_dn, FAM_ds, baseline_MG94xREV_omega, tau, value, gc_content, n_genes, ratio) %>% 
  mutate_all(~ as.numeric(.))

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
ds_tree <- dups_dnds %>% select(id, FAM_ds)
ds_tree <- merge(ds_tree, longest_transcript[c('YOgn', 'species')], by.x = 'id', by.y = 'YOgn')

tree <- ape::read.tree('./MEGA_Tree.nwk')

library(ggtree)
ggplot(ds_tree, aes(x = FAM_ds, y = factor(species, levels = rev(c('dmel','dyak','dana','dpse','dper','dvir','dmoj','dwil'))))) +
  geom_jitter() +
  xlim(2,0) +
  ylab('') +
  xlab('Ds') +
  theme_bw() +
  theme(panel.border = element_blank()) +
  ggtree(tree) + geom_tiplab(aes(label = label)) 

ggsave(filename = './Plots/ds_phylogeny.jpg', width = 13, height = 5)


# 
tree <- ape::read.tree('./MEGA_Tree.nwk')

library(ggtree)


node_data <- data.frame(
  node = tree['tip.label'],
  category1 = sample(1:3, length(tree['tip.label']), replace = TRUE),
  category2 = sample(1:3, length(tree['tip.label']), replace = TRUE)
)



func_counts <- all_func %>%
  mutate(spec = case_when(str_detect(id, "AN") ~ 'dana',
                          str_detect(id, "ME") ~ 'dmel',
                          str_detect(id, "MO") ~ 'dmoj',
                          str_detect(id, "PE") ~ 'dper',
                          str_detect(id, "PS") ~ 'dpse',
                          str_detect(id, "VI") ~ 'dvir',
                          str_detect(id, "WI") ~ 'dwil',
                          str_detect(id, "YA")~ 'dyak')) %>%
  select(spec, func_pseudo) %>%
  mutate(func_pseudo = gsub('_dup1', '', func_pseudo),
         func_pseudo = gsub('_dup2', '', func_pseudo),
         func_pseudo = gsub('neofunctionalized_dup1', 'neofunctionalized', func_pseudo),
         func_pseudo = gsub('neofunctionalized_dup2', 'neofunctionalized', func_pseudo),
         func_pseudo = gsub('pseudo_dup_1', 'pseudo', func_pseudo),
         func_pseudo = gsub('pseudo_dup_2', 'pseudo', func_pseudo),
         func_pseudo = gsub('sub', 'subfunctionalized', func_pseudo),
         func_pseudo = gsub('spec', 'specialized', func_pseudo),
         func_pseudo = gsub('cons', 'conserved', func_pseudo))

func_counts <- as.data.frame(table(func_counts$spec, func_counts$func_pseudo))
colnames(func_counts) <- c('spec','category','Freq')

# mech per species

all_mech <- dups_dnds %>%
  filter(!is.na(mech) & dup_or_ortho == 'dup')

mech_counts <- all_mech %>%
  mutate(spec = case_when(str_detect(id, "AN") ~ 'dana',
                          str_detect(id, "ME") ~ 'dmel',
                          str_detect(id, "MO") ~ 'dmoj',
                          str_detect(id, "PE") ~ 'dper',
                          str_detect(id, "PS") ~ 'dpse',
                          str_detect(id, "VI") ~ 'dvir',
                          str_detect(id, "WI") ~ 'dwil',
                          str_detect(id, "YA")~ 'dyak')) %>%
  select(spec,mech) %>%
  mutate(mech = gsub('_dup_1', '', mech),
         mech = gsub('_dup_2', '', mech))

mech_counts <- as.data.frame(table(mech_counts$spec, mech_counts$mech))
colnames(mech_counts) <- c('spec','category','Freq')




# plots 
ggtree_obj <- ggtree(tree) + geom_tiplab() 


func_levels <- c(rep('conserved', 8), rep('neofunctionalized', 8), 
                 rep('specialized', 8), rep('subfunctionalized', 8),
                 rep('pseudo', 8), rep('pseudo_both', 8))


func_bars <- ggplot(data = func_counts, aes(x = "", y = Freq, fill = category)) +
  geom_bar(width = 1, stat = "identity", position = 'fill') +
  geom_text(data = subset(func_counts, Freq != 0), aes(label = Freq), position = position_fill(vjust = 0.5), colour = 'black', size = 3) +  
  facet_grid(factor(spec, levels = c('dmel','dyak','dana','dpse','dper','dvir','dmoj','dwil')) ~ .) + 
  theme_void() +
  labs(title = '  Function') +
  scale_fill_manual(values = c('#6CBC4D','#E9AB17','#C05780','#778899','gray','darkgray'),
                    breaks = c('subfunctionalized','specialized','neofunctionalized','conserved','pseudo','pseudo_both')) +
  theme(legend.position = "bottom", legend.title = element_blank(), strip.text = element_blank()) +
  coord_flip()


mech_pies <- ggplot(data = mech_counts, aes(x = "", y = Freq, group = category, fill = category)) +
  geom_bar(width = 1, stat = "identity", position = position_fill()) +
  geom_text(aes(label = Freq), position = position_fill(vjust = 0.5), colour = 'black', size = 3) +
  coord_polar("y", start = 0) + 
  facet_grid(factor(spec, levels = c('dmel','dyak','dana','dpse','dper','dvir','dmoj','dwil')) ~ .) + 
  theme_void() +
  scale_fill_manual(values = c("#0079C2", "#FF828B", "808080")) + ###
  labs(title = 'Mechanism') +
  theme(legend.position = 'bottom', legend.title = element_blank(), strip.text = element_blank())

mech_func_phylogeny <- ggarrange(ggtree_obj, mech_pies, func_bars, 
                                 nrow = 1, ncol = 3, 
                                 widths = c(2, 0.39, 0.9),
                                 heights = c(1,1,1))

ggsave(mech_func_phylogeny, filename = './Plots/mech_func_phylogeny.jpg', width = 18, height = 7)


ggtree_obj2 <- ggtree_obj + xlim(0, 0.9)
ggsave(ggtree_obj2, filename = './Plots/tall_mega_phylogeny.jpg',width = 4, height = 7)


# sex bias
sex_bias <- read.csv("./sex_bias.tsv", sep = '')

orig_func <- orig_func %>%
  select(dup_1, dup_2, func = func_pseudo)

colnames(sex_bias) <- c('dup_1','dup_1_sex_bias')
sex_bias_dups <- merge(orig_func, sex_bias, by = 'dup_1')

colnames(sex_bias) <- c('dup_2','dup_2_sex_bias')
sex_bias_dups <- merge(sex_bias_dups, sex_bias, by = 'dup_2')


sex_bias <- as.data.frame(table(sex_bias_dups$func, sex_bias_dups$dup_1_sex_bias, sex_bias_dups$dup_2_sex_bias))


colnames(sex_bias) <- c('func','dup1_bias','dup2_bias','freq')




sex_bias <- sex_bias_dups %>%
  mutate(bias = case_when(
    dup_1_sex_bias == 'female_biased' & dup_2_sex_bias == 'female_biased' ~ 'both_fem',
    dup_1_sex_bias == 'male_biased' & dup_2_sex_bias == 'male_biased' ~ 'both_male',
    dup_1_sex_bias == 'neutral' & dup_2_sex_bias == 'neutral' ~ 'both_neut',
    
    dup_1_sex_bias == 'female_biased' & dup_2_sex_bias == 'male_biased' ~ 'switched_bias',
    dup_1_sex_bias == 'male_biased' & dup_2_sex_bias == 'female_biased' ~ 'switched_bias',

    dup_1_sex_bias == 'male_biased' & dup_2_sex_bias == 'neutral' ~ 'one_male',
    dup_1_sex_bias == 'neutral' & dup_2_sex_bias == 'male_biased' ~ 'one_male',
    
    dup_1_sex_bias == 'female_biased' & dup_2_sex_bias == 'neutral' ~ 'one_fem',
    dup_1_sex_bias == 'neutral' & dup_2_sex_bias == 'female_biased' ~ 'one_fem'
    
  ),
  
  sex_bias = case_when(
    bias == 'one_male' & func == 'neo_dup1' & dup_1_sex_bias == 'male_biased' ~ 'neo_new_M',
    bias == 'one_male' & func == 'neo_dup2' & dup_2_sex_bias == 'male_biased' ~ 'neo_new_M',
    
    bias == 'one_fem' & func == 'neo_dup1' & dup_1_sex_bias == 'female_biased' ~ 'neo_new_F',
    bias == 'one_fem' & func == 'neo_dup2' & dup_2_sex_bias == 'female_biased' ~ 'neo_new_F',
    
    
    bias == 'switched_bias' & func == 'neo_dup1' & dup_1_sex_bias == 'male_biased' ~ 'neo_switch_to_M_fromF',
    bias == 'switched_bias' & func == 'neo_dup2' & dup_2_sex_bias == 'male_biased' ~ 'neo_switch_to_M_fromF',
    bias == 'switched_bias' & func == 'neo_dup1' & dup_1_sex_bias == 'female_biased' ~ 'neo_switch_to_F_fromM',
    bias == 'switched_bias' & func == 'neo_dup2' & dup_2_sex_bias == 'female_biased' ~ 'neo_switch_to_F_fromM',
    bias == 'switched_bias' & func == 'specialized' ~ 'spec_switch_bias',
    T ~ bias)
  )

sex_bias <- sex_bias %>% mutate(func = gsub('neo_dup1', 'neo', func),
                  func = gsub('neo_dup2', 'neo', func))
sex_bias <- as.data.frame(table(sex_bias$sex_bias, sex_bias$func))
sex_bias <- pivot_wider(sex_bias, names_from = 'Var2', values_from = 'Freq')
sex_bias$Var1 <- c('both female','both male','both neutral','neo gained female bias',
            'neo gained male bias', 'neo switched from M to F', # no 'neo switched to F from M',
            'one female biased copy (the non-neo copy for neo)', 
            'one male biased copy (the non-neo one for neo)', 'switched bias')
colnames(sex_bias) <- c('bias','conserved','neo','specialized','sub')




# across species comparison

funcs <- read.csv("./CDROM_CLOUD_pseudo_funcs.tsv", sep = '')


funcs <- funcs %>%
  mutate(species = case_when(grepl('AN', dup_1) ~ 'dana',
                             grepl('ME', dup_1) ~ 'dmel',
                             grepl('MO', dup_1) ~ 'dmoj',
                             grepl('PE', dup_1) ~ 'dper',
                             grepl('PS', dup_1) ~ 'dpse',
                             grepl('VI', dup_1) ~ 'dvir',
                             grepl('WI', dup_1) ~ 'dwil',
                             grepl('YA', dup_1) ~ 'dyak')) %>%
  mutate(CLOUD_pseudo = case_when(str_detect(CLOUD_pseudo, 'pseudo') ~ 'pseudo', T ~ CLOUD_pseudo))  %>%
  mutate(CDROM_pseudo = case_when(str_detect(CDROM_pseudo, 'pseudo') ~ 'pseudo', T ~ CDROM_pseudo))  %>%
  mutate(CLOUD_pseudo = case_when(str_detect(CLOUD_pseudo, 'neo') ~ 'neo', T ~ CLOUD_pseudo))  %>%
  mutate(CDROM_pseudo = case_when(str_detect(CDROM_pseudo, 'neo') ~ 'neo', T ~ CDROM_pseudo))

funcs_counts <- as.data.frame(table(funcs$species, funcs$CLOUD_pseudo))
colnames(funcs_counts) <- c('spec','category','Freq')


CL <- funcs_counts %>% 
  filter(Freq != 0) %>%
  filter(category != 'pseudo') %>%
  ggplot(aes(x="", y=Freq, group=category, fill=category)) +
  geom_bar(width = 1, stat = "identity", position = position_fill()) +
  geom_text(aes(label = Freq), position = position_fill(vjust = 0.5), colour = 'black', size = 3) +
  coord_polar("y", start=0) + 
  labs(title = 'CLOUD') +
  facet_grid(.~factor(spec, levels= c('dmel','dyak','dana','dpse','dper','dvir','dmoj','dwil'))) + 
  theme_void() #+
#scale_fill_manual(values=c("#0079C2", "#FF828B", "808080")) + ###
#labs(title = 'Mechanism') +
#theme(legend.position = 'bottom', legend.title = element_blank(), strip.text = element_blank())




funcs_counts_CDROM <- as.data.frame(table(funcs$species, funcs$CDROM_pseudo))
colnames(funcs_counts_CDROM) <- c('spec','category','Freq')


CD <- funcs_counts_CDROM %>% 
  filter(Freq != 0) %>%
  filter(category != 'pseudo') %>%
  ggplot(aes(x="", y=Freq, group=category, fill=category)) +
  geom_bar(width = 1, stat = "identity", position = position_fill()) +
  geom_text(aes(label = Freq), position = position_fill(vjust = 0.5), colour = 'black', size = 3) +
  coord_polar("y", start=0) + 
  labs(title = 'CDROM') +
  facet_grid(.~factor(spec, levels= c('dmel','dyak','dana','dpse','dper','dvir','dmoj','dwil'))) + 
  theme_void() #+
#scale_fill_manual(values=c("#0079C2", "#FF828B", "808080")) + ###
#labs(title = 'Mechanism') +
#theme(legend.position = 'bottom', legend.title = element_blank(), strip.text = element_blank())






across_species <- ggarrange(CD, CL, nrow = 2, ncol = 1, align = 'hv')

ggsave(across_species, filename = './Plots/CDROM_CLOUD_across_species.jpg', width = 10, height = 4)




# compare CDROM and CLOUD classifications 
funcs_counts_both <- as.data.frame(table(funcs$species, funcs$CDROM_pseudo, funcs$CLOUD_pseudo))

colnames(funcs_counts_both) <- c('species','CDROM','CLOUD','count')


CDROM_v_CLOUD <- funcs_counts_both %>% 
  filter(CDROM != "pseudo" & CLOUD != "pseudo") %>%
  group_by(CDROM, CLOUD) %>%
  summarise(count = sum(count), .groups = 'drop') %>%
  
  ggplot(., aes(x = CDROM, y = CLOUD, fill = count)) +
  geom_tile(color = "white") +
  geom_text(aes(label = count), vjust = 1, color = "black") +
  scale_fill_gradient(low = "white", high = "steelblue") +
  theme_minimal()


ggsave(CDROM_v_CLOUD, filename = './Plots/CDROM_CLOUD_classification.jpg', width = 6, height = 5)



