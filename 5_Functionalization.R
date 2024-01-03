



source('./startup.R')


# read in duplicate genes with their ancestral copy 
dups <- read.csv("./Dup_Pairs_Ancestral.tsv", sep="")

# read in ortholog pairs 
ortho_pairs <- read.csv("./Ortholog_Pairs.tsv", sep="")

# read in expression data 
exp <- read.csv("./Expression_Data.tsv", sep="")




# separate each gene id (dup1,dup2,ancestral,single copy one, single copy two)

dup_1 <- dups[c('dup_1')]
dup_2 <- dups[c('dup_2')]
anc <- dups[c('ancestral_copy')]
sc_x <- ortho_pairs[c('YOgn.x')]
sc_y <- ortho_pairs[c('YOgn.y')]

# get expression data for each single copy gene
colnames(sc_x) <- 'YOgnID'
exp_sc_x <- merge(sc_x,exp,by='YOgnID')
colnames(sc_y) <- 'YOgnID'
exp_sc_y <- merge(sc_y,exp,by='YOgnID')
rm(sc_x, sc_y)

# get expression data for dups and anc
rownames(exp) <- exp$YOgnID
exp <- exp %>% select(-YOgnID)
exp_dup_1 <- exp[row.names(exp) %in% dup_1$dup_1, ]
exp_dup_2 <- exp[row.names(exp) %in% dup_2$dup_2, ]
exp_anc <- exp[row.names(exp) %in% anc$ancestral_copy, ]

# get combined dup expression
exp_dups_combined <- exp_dup_1 + exp_dup_2

# get relative expression levels 
exp_dup_1 <- exp_dup_1 / rowSums(exp_dup_1)
exp_dup_2 <- exp_dup_2 / rowSums(exp_dup_2)
exp_anc <- exp_anc / rowSums(exp_anc)
exp_dups_combined <- exp_dups_combined / rowSums(exp_dups_combined)

# get relative expression levels for single copy genes
rownames(exp_sc_x) <- make.names(exp_sc_x$YOgnID, unique=TRUE)
exp_sc_x <- exp_sc_x %>% select(-YOgnID)
rownames(exp_sc_y) <- make.names(exp_sc_y$YOgnID, unique=TRUE)
exp_sc_y <- exp_sc_y %>% select(-YOgnID)

exp_sc_x <- exp_sc_x / rowSums(exp_sc_x)
exp_sc_y <- exp_sc_y / rowSums(exp_sc_y)

# get euclidean distances (ed) values
dup_1_a <- (rowSums((exp_dup_1 - exp_anc) ^ 2)) ^ (1/2)
dup_2_a <- (rowSums((exp_dup_2 - exp_anc) ^ 2)) ^ (1/2)
dups_combined_a <- (rowSums((exp_dups_combined - exp_anc) ^ 2)) ^ (1/2)


## OLD script:







# get single copy ed values
sc$gn_2 <- make.names(sc$gn_2, unique=TRUE)
rownam <- rownames(exp_sc_y)
match_positions <- match(sc$gn_2, rownam)
rn <- rownam[match_positions]
exp_sc_y <- exp_sc_y[rn, ]

# calculate the sc euclidean distance
exp_sc <- exp_sc_x - exp_sc_y
exp_sc <- exp_sc ^ 2
exp_sc <- rowSums(exp_sc)
sc_ed <- exp_sc ^ (1/2)
sc_ed <- as.data.frame(sc_ed)

# classify into functional groups 
ed <- as.data.frame(cbind(dup_1_a, dup_2_a, dups_combined_a))

# remove duplicates with any non-expressed pairs
ed <- na.omit(ed)

# classify

sc_ed <- na.omit(sc_ed)




ed$conserv <- NA
ed$neo_dup1 <- NA
ed$neo_dup2 <- NA
ed$subfun <- NA
ed$specializ <- NA

sc_ed$sc_ed <- sc_ed$sc_ed / 2



for (row_num in 1:nrow(ed)){
  
  sample <- as.data.frame(sample(sc_ed$sc_ed,size=1000))
  colnames(sample) <- 'sc_ed'
  
  d1 <- ed[row_num,1]
  d2 <- ed[row_num,2]
  dd <- ed[row_num,3]
  
  ed[row_num,'conserv'] <- sum((d1 <= sample$sc_ed) & (d2 <= sample$sc_ed))
  ed[row_num,'neo_dup1'] <- sum((d1 > sample$sc_ed & d2 <= sample$sc_ed))
  ed[row_num,'neo_dup2'] <- sum((d1 <= sample$sc_ed & d2 > sample$sc_ed))
  ed[row_num,'subfun'] <- sum((d1 > sample$sc_ed & d2 > sample$sc_ed & dd <= sample$sc_ed))
  ed[row_num,'specializ'] <- sum((d1 > sample$sc_ed & d2 > sample$sc_ed & dd > sample$sc_ed))
  
}


ed$func <- names(ed[, 4:8])[apply(ed[, 4:8], 1, which.max)]

table(ed$func)

# not divided by 2:
#  conserv  neo_dup1  neo_dup2 specializ    subfun 
# 266       158       144       610         3 


# format results
ed$dup_1 <- rownames(ed)
ed <- merge(ed,dups,by='dup_1')
ed <- ed[,c(1,11,12,2:4,10,5:8)]

# write to file








