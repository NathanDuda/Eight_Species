


source('./startup.R')

# import all annotations from YO 
all_annotations <- read.delim("C:/Users/17735/Downloads/Eight_Species/Raw_YO_Data/YO_Annotations/all.YO.gtf", header=FALSE)
colnames(all_annotations) <- c('chrom','StringTie_or_FlyBase','type','start','end','x','strand','y','id')

# format the transcripts
transcripts <- all_annotations[all_annotations$type=='transcript',]
transcripts <- transcripts %>% separate(id, into = c('YOtr','YOgn','old','tss'), sep = ";")

# remove the YOtrs that are inside other YOtrs because won't be the longest
transcripts <- transcripts[!grepl('contained_in YOtr',transcripts$tss),]

# keep only the longest YOtr per YOgn
transcripts$length <- transcripts$end - transcripts$start

transcripts <- transcripts %>% 
  group_by(YOgn) %>%
  filter(length == max(length))
longest_transcript <- transcripts[!duplicated(transcripts[c('YOgn')]),]


# notes:
## MSTRG is the default gene name from StringTie
## if StringTie, its a new annotation from YO
## tss_id is present when the transcript starts at the same location as the first exon


# format the longest_transcript dataframe 
longest_transcript$YOtr <- gsub('transcript_id ','',longest_transcript$YOtr)
longest_transcript$YOgn <- gsub(' gene_id ','',longest_transcript$YOgn)
longest_transcript$old <- gsub(' oId ','',longest_transcript$old)

longest_transcript[grepl('MSTRG', longest_transcript$old), 'old'] <- NA

longest_transcript <- longest_transcript[c('YOgn','YOtr','old','chrom','strand','start','end')]


# get the sequences for 




