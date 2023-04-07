# CRISPR:

crispr <- read.table("Data/CRISPR/crispr_pilecr.tsv", header=TRUE, sep="\t")
mag_info <- read.table("Data/FinalBinSet/MAG_Characteristics.tsv", header=TRUE, sep="\t")

library(tidyverse)

crispr <- left_join(crispr, mag_info, by=c("Genome"="genome"))
crispr <- crispr %>% filter(Putative_CRISPR_arrays_pilecr > 0)

ggplot(crispr %>% arrange(desc(Putative_CRISPR_arrays_pilecr)), aes(x=Genus, y=Putative_CRISPR_arrays_pilecr))+
  geom_col()+
  facet_grid(.~Phylum, scales="free")+
  theme_bw()+
  #my_fave_theme+
  theme(strip.background = element_blank(),
        strip.text = element_text(face="bold", color="black", size=4),
        axis.text.x = element_text(angle=90))

ggsave("Figures/CRISPR.pdf", width = 11, height=5)
