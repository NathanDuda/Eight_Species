

library(tidyverse)
library(ape)



rel_exp <- read.csv("./Relative_Expression.tsv", sep="")
colnames(rel_exp)[1] <- 'id'


#devtools::install_github("hr1912/TreeExp")

orthogroups <- read.delim2("./OrthoFinder_Output/Results_Jan01/Orthogroups/Orthogroups.tsv")


orthogroups <- orthogroups %>%
  pivot_longer(cols = c(2:ncol(orthogroups))) %>%
  separate_rows(value, sep = ', ') %>%
  mutate_all(~ifelse(. == "", NA, .))



og <- 'OG0000968'

for (og in unique(orthogroups$Orthogroup)) {
  
  og_gene_exp <- orthogroups %>%
    filter(Orthogroup == og) %>% # keep only genes for given orthogroup 
    na.omit() %>%
    select(id = value) %>%
    left_join(., rel_exp, by = 'id') # merge with expression values 
  
  exp_cartesian_product <- merge(og_gene_exp, og_gene_exp, by = NULL)
  
  # calculate euclidean distances for each cartesian product pair (get euclidean distance matrix)
  exp.y <- exp_cartesian_product %>%
    select(ends_with('.y')) %>%
    mutate(id.y = make.unique(id.y)) %>%
    column_to_rownames('id.y')
  
  exp.x <- exp_cartesian_product %>%
    select(ends_with('.x')) %>%
    mutate(id.x = make.unique(id.x)) %>%
    column_to_rownames('id.x')
  
  euc_dist <- (rowSums((exp.x - exp.y) ^ 2)) ^ (1/2) # euclidean distance formula
  
  #
  ed_matrix <- exp_cartesian_product %>%
    select(id.x, id.y) %>%
    mutate(ed = euc_dist) %>%
    filter(id.x != id.y) %>% # remove eds between the same gene 
    mutate(id_min = pmin(id.x, id.y), id_max = pmax(id.x, id.y)) %>% # remove reciprocally identical eds
    distinct(id_min, id_max, .keep_all = TRUE) %>%
    select(-id_min, -id_max)
  
  
  
  expression_matrix <- ed_matrix %>%
    pivot_wider(names_from = 'id.y', values_from = 'ed') %>%
    column_to_rownames('id.x') %>%
    add_row(., .before = 1)
  rownames(expression_matrix)[1:nrow(expression_matrix) - 1] <- colnames(expression_matrix)
  expression_matrix$new <- NA
  colnames(expression_matrix) <- rownames(expression_matrix)
  
  

  expression_tree <- midpoint(nj(as.dist(expression_matrix))) # midpoint-rooted neighbor-joining tree
  
  
  
  gene_tree <- read.tree(file = paste0('./OrthoFinder_Output/Results_Jan01/Resolved_Gene_Trees/', og, '_tree.txt'))
  gene_tree$tip.label <- sub(".*YOgn", "YOgn", gene_tree$tip.label) # remove the characters before the gene id 
  
  
  plot(expression_tree)
  plot(gene_tree)
  
  
  gene_matrix <- cophenetic(gene_tree)
  gene_matrix[upper.tri(gene_matrix, diag = T)] <- NA # keep only the lower half of the gene matrix 
  
  
  
  
  expression_matrix <- as.data.frame(expression_matrix)
  gene_matrix <- as.data.frame(gene_matrix)
  
  
  comparePhylo(gene_tree, expression_tree, plot = T)
  
  
  all.equal(gene_tree, expression_tree)
  
  treedist(gene_tree, expression_tree)
  
  
  
  
  
  gn_matrix <- gene_matrix %>%
    rownames_to_column('id.y') %>%
    pivot_longer(cols = 2:ncol(.)) %>%
    na.omit() %>%
    select(id.y, id.x = name, ed = value) %>%
    mutate(type = 'gene')
  
  ed_matrix <- ed_matrix %>%
    mutate(type = 'expression')
  
  t <- rbind(ed_matrix, gn_matrix)
  
  ggplot(t, aes(x = id.y, y = id.x, fill = ed)) +
    geom_tile() +
    facet_grid(.~type)
  
  

  t %>%
    mutate(pair = rep(1:(nrow(.) / 2), 2)) %>%
    ggplot(aes(x = type, y = ed)) +
      geom_point() +
      geom_line(aes(group=pair)) 
  
  
  
  
  
}










