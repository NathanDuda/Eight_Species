# Author: Nathan Duda
# Purpose: 
#   This script organizes the raw annotations from Yang & Oliver (YO).
#   Nucleotide and protein sequences are extracted for these annotations from
#   the FlyBase_2017_03_Genome.


source('./startup.R')

# import all annotations from YO 
all_annotations <- read.delim("./Raw_Data/YO_Annotations/all.YO.gtf", header=FALSE)
colnames(all_annotations) <- c('chrom','StringTie_or_FlyBase','type','start','end','x','strand','y','id')

# format the transcripts
transcripts <- all_annotations %>%
  filter(type == 'transcript') %>%
  separate(id, into = c('YOtr','YOgn','old','tss'), sep = "; ")

# remove the YOtrs that are inside other YOtrs because won't be the longest
transcripts <- transcripts %>%
  filter(!grepl('contained_in YOtr', tss))

# keep only the longest YOtr per YOgn (can be from StringTie or FlyBase)
longest_transcript <- transcripts %>%
  mutate(length = end - start) %>% 
  group_by(YOgn) %>%
  filter(length == max(length)) %>%
  distinct(YOgn, .keep_all = TRUE)

# notes:
## MSTRG is the default gene name from StringTie
## if StringTie, its a new annotation from YO
## tss_id is present when the transcript starts at the same location as the first exon

# format the longest_transcript dataframe 
longest_transcript <- longest_transcript %>%
  mutate(YOtr = gsub('transcript_id ', '', YOtr),
         YOgn = gsub(' gene_id ', '', YOgn),
         old = gsub(' oId ','', old),
         old = ifelse(grepl('MSTRG', old), NA, old)) %>%
  select(YOgn, YOtr, old, chrom, strand, start, end)


# get the sequences for the transcripts
fasta_sequences <- readDNAStringSet('./Raw_Data/FlyBase_2017_03_Genomes/all_chromosome.fasta')

species <- sub(".*species=([^;]+).*", "\\1", names(fasta_sequences))
scaffolds <- gsub(" .*", "", names(fasta_sequences))
names(fasta_sequences) <- paste0(species,'_',scaffolds)

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
    if (strand == "-") {
      seq <- reverseComplement(seq)}
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
                                                         longest_orf$start, 
                                                         longest_orf$end)
  }
  
  print(row_num)
}


# remove the ORFs with 'N' in the sequences and those without ORFs 
longest_transcript <- longest_transcript %>%
  filter(!grepl('N',longest_ORF),
         !is.na(longest_ORF)) %>%
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
  species_data <- longest_transcript %>%
    filter(species == species_name)
  
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


