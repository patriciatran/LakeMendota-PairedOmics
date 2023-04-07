samples.metadata <- read.csv("Data/Samples-mendota.csv")
jgi.ids <- read.csv("Data/SampleChar.csv")
jgi.ids$Lake <- "Lake Mendota"


library(tidyverse)

tests <- samples.metadata %>% 
  gather(Sample.type, Value, Metagenome:Metatranscriptome) %>%
  mutate(Sample.type = ifelse(Sample.type =="Metagenome",yes="Total metagenome",Sample.type)) %>%
  full_join(jgi.ids) %>% 
  dplyr::select(-Sample.collection.date) %>%
  arrange(Date.YYYY.MM.DD, Sample.Depth.meters)

write.table(tests, "Data/SupplementaryTable1.tsv", sep="\t", quote=FALSE, row.names = FALSE)
