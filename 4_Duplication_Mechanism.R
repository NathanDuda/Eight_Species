


source('./startup.R')


# read in duplicate genes with their ancestral copy 
dups <- read.csv("./Dup_Pairs_Ancestral.tsv", sep="")

# read in and format all gene annotations with exons 
all_annotations <- read.delim("./Raw_Data/YO_Annotations/all.YO.gtf", header=FALSE)
colnames(all_annotations) <- c('chrom','StringTie_or_FlyBase','type','start','end','x','strand','y','id')
all_annotations <- all_annotations %>% separate(id, into = c('YOtr','YOgn','old','tss'), sep = ";")

exon_counts <- all_annotations %>%
  # get unique exons to count 
  filter(type == 'exon') %>%
  distinct(chrom,start,end, .keep_all = T) %>%

  # count exons per YOtr
  group_by(YOtr) %>%
  mutate(n_exons = n()) %>%

  # format into only YOtr and exon amount
  ungroup() %>%
  select(YOtr,n_exons) %>%
  mutate(YOtr = gsub('transcript_id ', '', YOtr)) %>%
  distinct()

# get the YOgns for the YOtrs used 
longest_transcript <- read.csv("./longest_transcript.tsv", sep="")
yogn_yotr <- longest_transcript[c('YOgn','YOtr')]

# get the exon counts corresponding to the longest transcript used
exon_counts <- merge(yogn_yotr,exon_counts,by='YOtr')
exon_counts <- exon_counts[c('YOgn','n_exons')]

# merge exon counts with duplicate genes 
colnames(exon_counts) <- c('dup_1','dup_1_n_exons')
dup_exons <- left_join(dups,exon_counts,by='dup_1')

colnames(exon_counts) <- c('dup_2','dup_2_n_exons')
dup_exons <- left_join(dup_exons,exon_counts,by='dup_2')

colnames(exon_counts) <- c('ancestral_copy','ancestral_copy_n_exons')
dup_exons <- left_join(dup_exons,exon_counts,by='ancestral_copy')

# NA when diff transcript (not longest one) for same gene has same location 
# change these to 1 
dup_exons <- dup_exons %>% mutate_all(~ ifelse(is.na(.), 1, .))

# classify as DNA- or RNA-mediated based on exon number
dup_exons <- dup_exons %>%
  mutate(mech = case_when(dup_1_n_exons > 1 & dup_2_n_exons == 1 ~ 'RNA',
                          dup_2_n_exons > 1 & dup_1_n_exons == 1 ~ 'RNA',
                          dup_1_n_exons > 1 & dup_2_n_exons > 1 ~ 'DNA',
                          dup_2_n_exons > 1 & dup_1_n_exons > 1 ~ 'DNA',
                          dup_2_n_exons == 1 & dup_1_n_exons == 1 ~ 'unknown'))

# write the duplication mechanisms of each duplicate pair to file
write.table(dup_exons,file= 'Dup_Mechanism.tsv')

