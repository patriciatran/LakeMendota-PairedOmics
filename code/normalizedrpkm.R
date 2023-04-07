# Normalizing the RPKM data

library(tidyverse)

MAG.abundance <- read.table("Data/MAGS-abundance-long-format.tsv", header=TRUE, sep="\t")

MAG.abundance.matrix <- MAG.abundance %>% select(Genome, SAMPLE.NAME, Coverage) %>% spread("SAMPLE.NAME","Coverage") 

write.table(MAG.abundance.matrix, "Data/MAG_abundance_matrix_rel_abund.tsv", sep="\t", row.names = FALSE)

# RPKM normalized table:
phages.rpkm <- read.table("Data/phages_abundance_matrix_metatranscriptomes_rpkm.tsv", sep="\t", header=TRUE)
IS_reads <- read.table("Data/mt_read_counts_IS-1.tsv", sep="\t", header=TRUE)


phages.rpkm.gather <- phages.rpkm %>% gather("mt.name", "rpkm", 2:17)
phages.rpkm.gather$mt.name <- str_replace(phages.rpkm.gather$mt.name,"X","")
phages.rpkm.gather$mt.name <- str_replace(phages.rpkm.gather$mt.name,".07.","-07-")
phages.rpkm.gather$mt.name <- str_replace(phages.rpkm.gather$mt.name,".08.","-08-")
phages.rpkm.gather$mt.name <- str_replace(phages.rpkm.gather$mt.name,".09.","-09-")
phages.rpkm.gather$mt.name <- str_replace(phages.rpkm.gather$mt.name,".10.","-10-")
phages.rpkm.gather$mt.name <- str_replace(phages.rpkm.gather$mt.name,"-05-10-","-05_10m")
phages.rpkm.gather$mt.name <- str_replace(phages.rpkm.gather$mt.name,"25-10-","25_10m")

unique(phages.rpkm.gather$mt.name)
unique(IS_reads$mt.name)

phages.rpkm.gather <- left_join(phages.rpkm.gather, IS_reads %>% select(mt.name,IS_reads))

phages.rpkm.gather.norm <- phages.rpkm.gather %>% mutate(rpkm_normalized = rpkm/IS_reads)

phages.rpkm.gather.norm.matrix <- phages.rpkm.gather.norm %>% select(Genome, mt.name, rpkm_normalized) %>%
  spread(mt.name, rpkm_normalized)

write.table(phages.rpkm.gather.norm.matrix, "Data/phages_rpkm_normalized_matrix.tsv",
            sep="\t",quote = FALSE, row.names = FALSE)
