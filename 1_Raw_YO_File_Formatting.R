


source('./startup.R')

# import all annotations from YO 
all_annotations <- read.delim("C:/Users/17735/Downloads/Eight_Species/Raw_Data/YO_Annotations/all.YO.gtf", header=FALSE)
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


# get the sequences for the transcripts

fasta_sequences <- readDNAStringSet('./Raw_Data/FlyBase_2017_03_Genomes/all_chromosome.fasta')

species <- sub(".*species=([^;]+).*", "\\1", names(fasta_sequences))
scaffolds <- gsub(" .*", "", names(fasta_sequences))

names(fasta_sequences) <- paste0(species,'_',scaffolds)


longest_transcript$chrom <- longest_transcript$chrom

longest_transcript <- longest_transcript %>%
  mutate(chrom = case_when(grepl('AN',YOgn) ~ paste0('Dana_',chrom),
                           grepl('ME',YOgn) ~ paste0('Dmel_',chrom),
                           grepl('MO',YOgn) ~ paste0('Dmoj_',chrom),
                           grepl('PE',YOgn) ~ paste0('Dper_',chrom),
                           grepl('PS',YOgn) ~ paste0('Dpse_',chrom),
                           grepl('VI',YOgn) ~ paste0('Dvir_',chrom),
                           grepl('WI',YOgn) ~ paste0('Dwil_',chrom),
                           grepl('YA',YOgn) ~ paste0('Dyak_',chrom)))


# function to extract sequences
extract_sequence <- function(chrom, start, end, strand, fasta_sequences) {
  matches <- fasta_sequences[names(fasta_sequences) == chrom]
  if (length(matches) > 0) {
    seq <- subseq(matches, start, end)
    if (strand == "-") {seq <- reverseComplement(seq)}
    as.character(seq)
  } 
  else {NA}
}


# apply the function to create column 'extracted_sequence'
longest_transcript$extracted_sequence <- mapply(
  extract_sequence,
  longest_transcript$chrom,
  longest_transcript$start,
  longest_transcript$end,
  longest_transcript$strand,
  MoreArgs = list(fasta_sequences)
)


# get the longest ORFs 

longest_transcript$longest_ORF <- NA
library(ORFik)
for (row_num in 1:nrow(longest_transcript)) {
  row <- longest_transcript[row_num,]
  
  # find the ORFs in the given transcript 
  orfs <- as.data.frame(findORFs(row$extracted_sequence, start = "ATG"))
  
  if (nrow(orfs) > 0) {
    # keep only longest ORF and if multiple have same length, keep only first one
    longest_orf <- orfs[which.max(orfs$width), c("start", "end"), drop = FALSE]
  
    # place the ORF sequence into the longest_transcript dataframe 
    longest_transcript[row_num, 'longest_ORF'] <- substr(longest_transcript[row_num,'extracted_sequence'], 
                                                         longest_orf$start, longest_orf$end)
  }
  
  print(row_num)
}


# remove the ORFs with 'N' in the sequences and those without ORFs 
longest_transcript <- longest_transcript %>%
  filter(!grepl('N',longest_ORF)) %>% filter(!is.na(longest_ORF)) %>%
  rowwise() %>%
  mutate(prot = as.character(Biostrings::translate(DNAString(longest_ORF))))

# keep only sequences with protein length >=10
longest_transcript <- longest_transcript %>%
  mutate(prot_length = nchar(prot)) %>%
  filter(prot_length >= 10)

# add a species column 
longest_transcript <- longest_transcript %>%
  mutate(species = case_when(grepl('AN',YOgn) ~ 'dana',
                             grepl('ME',YOgn) ~ 'dmel',
                             grepl('MO',YOgn) ~ 'dmoj',
                             grepl('PE',YOgn) ~ 'dper',
                             grepl('PS',YOgn) ~ 'dpse',
                             grepl('VI',YOgn) ~ 'dvir',
                             grepl('WI',YOgn) ~ 'dwil',
                             grepl('YA',YOgn) ~ 'dyak'))

# write table to file 
write.table(longest_transcript,file='longest_transcript.tsv')



# write fasta files 
unique_species <- unique(longest_transcript$species)

for (species_name in unique_species) {
  species_data <- longest_transcript[longest_transcript$species == species_name, ]
  
  # for proteins:
  sequences <- paste0(">", species_data$YOgn, "\n", species_data$prot)
  fasta_content <- paste(sequences, collapse = "\n")
  
  file_path <- file.path('./Protein_Fastas/', paste0(species_name, "_prot.fasta"))
  writeLines(fasta_content, con = file_path)
  
  # for nucleotides:
  sequences <- paste0(">", species_data$YOgn, "\n", species_data$longest_ORF)
  fasta_content <- paste(sequences, collapse = "\n")
  
  file_path <- file.path('./Nucleotide_Fastas/', paste0(species_name, "_nuc.fasta"))
  writeLines(fasta_content, con = file_path)
}


