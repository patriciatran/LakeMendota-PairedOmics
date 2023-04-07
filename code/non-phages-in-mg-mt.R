# comparison of sourmash results of non-phages
# pulled out non phages sequences from the metagenoms and viromes and then ran SOURMASH on the genbank file

sourmash.res.mg <- read.csv("Data/Non_phages_in_metagenomes_against_genbank.csv")
sourmash.res.virome <- read.csv("Data/Non_phages_against_genbank.csv")


library(tidyverse)

sourmash.non.phages <- full_join(sourmash.res.mg, sourmash.res.virome)

sourmash.non.phages %>% 
  mutate(Phylum_Class_Order = paste(phylum,"_",class,"_", order)) %>%
  group_by(Phylum_Class_Order) %>% 
  mutate(Grouped_count = sum(count)) %>% 
  select(superkingdom,Phylum_Class_Order, Grouped_count, filename) %>% distinct() %>% arrange(desc(Grouped_count)) %>%
  ggplot(aes(x=Phylum_Class_Order, y=Grouped_count, fill=filename))+
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        legend.position = "bottom",
        strip.background = element_blank())+
  scale_fill_manual(values=c("turquoise","blue"))+
  facet_grid(.~superkingdom, scales="free", space = "free" )

ggsave("Figures/What_are_the_non_phages_2.pdf", width = 11, height = 8.5)

sourmash.compared <- sourmash.non.phages %>% mutate(Phylum_Class_Order = paste(phylum,"_",class,"_", order)) %>%
  group_by(Phylum_Class_Order) %>% 
  mutate(Grouped_count = sum(count)) %>% 
  select(superkingdom,Phylum_Class_Order, Grouped_count, filename) %>%
  distinct() %>% 
  arrange(desc(Grouped_count))

sourmash.compared

test1 <- sourmash.non.phages %>% select(filename, superkingdom, phylum, count) %>%
  group_by(filename, superkingdom, phylum) %>% 
  mutate(sumCount = sum(count))


test1

mg.sourmash <- test1 %>% select(-count)%>% distinct() %>%
  filter(filename=="Non_phages_in_metagenomes.sourmash.sketch") %>%
  mutate(new_name = paste0(superkingdom, phylum, sep=" "))

vir.sourmash <- test1 %>% select(-count)%>%
  distinct() %>%
  filter(filename=="Non_phages_in_viromes.sourmash.sketch")  %>%
  mutate(new_name = paste0(superkingdom, phylum, sep=" "))

pdf("Figures/NonPhages.in.samples.pdf", width = 10, height = 6)
par(mfrow=c(1,2))
pie(mg.sourmash$sumCount, labels=mg.sourmash$new_name, main="Non-phages in metagenomes")
pie(vir.sourmash$sumCount, labels=vir.sourmash$new_name, main="Non-phages in viromes")
dev.off()
