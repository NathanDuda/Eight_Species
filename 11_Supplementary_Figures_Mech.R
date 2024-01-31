


source('./startup.R')


mech <- read.csv("./Dup_Mechanism.tsv", sep="")

table(mech$mech)
table(mech$mech_anc_more_than_one)
table(mech$mech_anc_same_as_parental)


# regular mech pie 
mech_pie <- mech[c('mech')]

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

ggsave(filename = './Supplementary_Figures/regular_mech_pie.jpg', width = 5, height = 3, device='tiff', dpi=300)




# semiconservative mech pie 
mech_pie <- mech[c('mech_anc_more_than_one')]

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

ggsave(filename = './Supplementary_Figures/semiconservative_mech_pie.jpg', width = 5, height = 3, device='tiff', dpi=300)



# ultraconservative mech pie 
mech_pie <- mech[c('mech_anc_same_as_parental')]

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

ggsave(filename = './Supplementary_Figures/ultraconservative_mech_pie.jpg', width = 5, height = 3, device='tiff', dpi=300)






