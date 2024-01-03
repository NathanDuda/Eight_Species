
# this script also gets all ortholog pairs of each duplicate pair 


source('./startup.R')

# read in duplicate pairs 
dups <- read.csv("./Duplicate_Pairs.tsv", sep="")

# read in orthologs of duplicate pairs 
orthologs <- read.csv2("./Dup_Pair_Orthologs.tsv", sep="")

# read in expression data 
expression <- read.csv("./Expression_Data.tsv", sep="")

# clean the column names of the orthologs dataframe and replace empty cells with NAs
orthologs <- orthologs %>% 
  rename_all(~gsub("_prot", "", .)) %>%
  mutate_all(~na_if(.,""))

# define the closest species of each species using a dictionary 
library(hash)
closest_species_dict <- hash(
  dana = 'dpse',
  dmel = 'dyak',
  dmoj = 'dvir',
  dper = 'dpse',
  dpse = 'dper',
  dvir = 'dmoj',
  dwil = 'dvir',
  dyak = 'dmel'
)

# create a list of species of how related they are to each other going down the phylogeny 
ordered_species <- c('dyak','dmel','dana','dpse','dper','dwil','dvir','dmoj',
                     'dvir','dwil','dper','dpse','dana','dmel','dyak')

# create function to find the closest expressed ortholog to each duplicate pair
find_closest_ortholog <- function(row,species) {
  closest_species <- closest_species_dict[[species]]
  closest_gene <- row[[closest_species]]
  
  if (is.na(closest_gene) | grepl(',',closest_gene) | (!closest_gene %in% expression$YOgnID)){
    n = 0
    non_expressed = 0 
    
    index <- min(which(ordered_species == closest_species))
    
    while(is.na(closest_gene) | grepl(',',closest_gene) | (!closest_gene %in% expression$YOgnID)){
      n = n + 1
      non_expressed = non_expressed + 1
      if ((index + n) > (length(ordered_species))) {return(NA)}
      closest_species <- ordered_species[index+n]
      closest_gene <- row[[closest_species]]
    }
  } 
  return(closest_gene)
}


# apply the function to each row of the ortholog table 
orthologs$ancestral_copy <- NA
for (row_num in 1:nrow(orthologs)) {
  row <- orthologs[row_num,] 
  orthologs[row_num,'ancestral_copy'] <- find_closest_ortholog(row, species = row$duplicate_pair_species)
}

# merge the ancestral copy with the duplicate pairs
ancestral_copy <- orthologs[c('Orthogroup','ancestral_copy')]
dups <- merge(dups,ancestral_copy,by='Orthogroup')

# remove duplicates without ancestral copies (the orthologs they had weren't expressed)
dups <- na.omit(dups)

# write the duplicate pairs with their ancestral copy to file
write.table(dups,'Dup_Pairs_Ancestral.tsv')



### all combinations of ortholog pairs

all_ortholog_pairs <- orthologs[c(1:9)] %>%
  mutate_all(~ifelse(grepl(",", .), NA, .)) %>%
  pivot_longer(cols = c(2:9))

all_ortholog_pairs <- 
  merge(all_ortholog_pairs,all_ortholog_pairs,by='Orthogroup') %>%
  na.omit() %>%
  filter(name.x!=name.y) %>% 

  # keep only one of the reciprocal combination pairs  
  mutate(x_y = paste0(value.x,'_',value.y)) %>%
  mutate(y_x = paste0(value.y,'_',value.x)) %>%
  filter(!(x_y > y_x)) %>%
  select(-x_y,-y_x)


colnames(all_ortholog_pairs) <- c('Orthogroup','species.x','YOgn.x','species.y','YOgn.y')


# keep orthologs with expression data and are expressed in at least one tissue
expressed_ortholog_pairs <- all_ortholog_pairs %>%
  filter((YOgn.x %in% expression$YOgnID) & (YOgn.y %in% expression$YOgnID))

# write ortholog pairs to file
write.table(expressed_ortholog_pairs,file = './Ortholog_Pairs.tsv')


