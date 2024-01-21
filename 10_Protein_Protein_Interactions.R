


source('./startup.R')


# import the three ppi data
ppi_1 <- read.delim("./PPI_Data/FlyBi_data_20210205.csv")
ppi_1 <- ppi_1[c('FBgn1','FBgn2')]

ppi_2 <- read.delim("./PPI_Data/fly_other_physical.txt")
ppi_2 <- ppi_2[c('FLY_GENE1','FLY_GENE2')]

ppi_3 <- read.delim("./PPI_Data/flybase_ppi.txt")
ppi_3 <- ppi_3[c('FLY_GENE1','FLY_GENE2')]

# combine all three datasets
ppi_2 <- rbind(ppi_2,ppi_3)
colnames(ppi_2) <- colnames(ppi_1)
ppi <- rbind(ppi_1,ppi_2)


ppi <- ppi[!duplicated(ppi),]

ppi_1 <- ppi[1]
ppi_2 <- ppi[2]

colnames(ppi_1) <- 'ppi'
colnames(ppi_2) <- 'ppi'

ppi <- rbind(ppi_1,ppi_2)

ppi <- as.data.frame(table(ppi$ppi))




fbgn_fbtr <- read.delim("C:/Users/17735/Downloads/Eight_Species/PPI_Data/dmel-all-r6.16.gtf", header=FALSE)


fbgn_fbtr$FBgn <- str_extract(fbgn_fbtr$V9, "(?<=gene_id\\s).*?(?=;\\s)")
fbgn_fbtr$FBtr <- str_extract(fbgn_fbtr$V9, "(?<=transcript_id\\s).*?(?=;\\s)")

fbgn_fbtr <- fbgn_fbtr[c('FBgn','FBtr')]

fbgn_fbtr <- na.omit(fbgn_fbtr)
fbgn_fbtr <- fbgn_fbtr[!duplicated(fbgn_fbtr),]




longest_transcript <- read.csv("./longest_transcript.tsv", sep="")
longest_transcript <- longest_transcript[longest_transcript$species=='dmel',]

longest_transcript <- longest_transcript %>%
  mutate(FBtr = case_when(grepl('FBtr',old) ~ old,
                          grepl('FBgn',old) ~ NA)) %>%
  left_join(.,fbgn_fbtr,by='FBtr') %>%
  mutate(FBgn = case_when(grepl('FBgn',old) ~ old,
                          grepl('FBtr',old) ~ FBgn))

yogn_ppi <- merge(longest_transcript,ppi,by.x='FBgn',by.y='Var1')
yogn_ppi <- yogn_ppi[c('YOgn','Freq')]
colnames(yogn_ppi) <- c('YOgnID','ppi')


all_plotting_categories <- read.csv("./all_plotting_categories.tsv", sep="")
all_plotting_categories <- all_plotting_categories[grepl('ME',all_plotting_categories$YOgnID),]




ppi_plot <- left_join(all_plotting_categories,yogn_ppi,by='YOgnID')

#ppi_plot$ppi[is.na(ppi_plot$ppi)] <- 0


func_levels <- c('All_Ortho','All_Dup','DNA','RNA(new)','RNA(anc)','Cons','Neo(new)','Neo(anc)','Spec','Sub')

ggplot(ppi_plot, aes(x=factor(category, levels = func_levels),y=ppi)) +
  geom_boxplot(color=c('#778899','#D4E056','#F67280',"#E9AB17","#E9AB17","#3B7BBD",'#E23F51','#E23F51','#6CBC4D','#F18244')) + 
  ylab('PPI') +
  xlab('') +
  ylim(0,150) +
  stat_compare_means(aes(label = after_stat(p.signif)),method = "t.test", 
                     ref.group = 'All_Ortho',hide.ns = TRUE, size=5) +
  geom_hline(yintercept=median(na.omit((ppi_plot[ppi_plot$category=='All_Ortho',])['ppi']$ppi)),
             linetype="dashed", color = '#778899',linewidth=0.5) +
  theme_bw()

ggsave(filename = './Plots/ppi_boxplot.jpg', width = 7, height = 5)








