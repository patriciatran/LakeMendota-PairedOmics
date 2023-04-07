# Virus Taxonomy:
library(tidyverse)

#Load samples
virus.taxo <- read.table("Data/Cat_phages_all_metaG_viromes.VIVID.virus-taxonomy.tsv", header=TRUE, sep="\t")
sample.match <- read.csv("Data/GaID_match.csv")

head(virus.taxo)

# Match

virus.taxo$GaID <- str_replace(virus.taxo$scaffold, "_.*","")

virus.taxo <- left_join(virus.taxo, sample.match)
library(lubridate)


# Plot by Order:
virus.taxo %>% group_by(Date, order, TYPE, Depth) %>% tally() %>%
  ggplot(aes(x=mdy(Date), y=n, fill=order))+
  geom_col()+
  facet_grid(Depth~TYPE)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())

#  Plot by family:
virus.taxo %>% group_by(Date, order, family, TYPE, Depth) %>% tally() %>%
  mutate(new_family = ifelse(grepl("Myoriridae|Siphoviridae|Podoviridae|unknown|ambiguous|unassigned", family),
                             family,
                             "Other viral taxa")) %>%
  mutate(new_family_order = ifelse(grepl("Other", new_family),
                                   "Other taxa",
                                   paste0(order,"(",family,")"))) %>%
  ggplot(aes(x=mdy(Date), y=n, fill=new_family_order))+
  geom_bar(color="black", stat="identity", position="fill")+
  #geom_col(color="black")+
  facet_grid(Depth~TYPE)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())+
  ylab("Phages count")+
  xlab("Sample date")+
  scale_fill_manual(values=c("#49A5D3",
                             "#EE461E","#EE921E","#FFEC1D","#AC6D21",
                             "#8B6EB5",
                             "#6EB57D",
                             "dark grey","dark grey",
                             "white"))
#  scale_fill_manual(values=c("red","orange","yellow","brown","pink","grey","white"),
#                    name="Viral taxonomy \nOrder and family")

ggsave("Figures/Viral_taxonomy__Order_family.pdf", width = 11, height = 8.5)

virus.taxo %>% group_by(order) %>% tally() %>% arrange(desc(n))

# Ok there aren't a lot of difference so let's see by genomem quality:
vibrant.qual.mg <- read.table("Data/VIBRANT_genome_quality_Cat_Assemblies_Metagenomes.tsv", sep="\t", header=TRUE)
vibrant.qual.virome <- read.table("Data/VIBRANT_genome_quality_Cat_Assemblies_Viromes.tsv", sep="\t", header=TRUE)


#combine files
vibrant.quality.all <- rbind(vibrant.qual.mg, vibrant.qual.virome)



virus.taxo <- left_join(virus.taxo, vibrant.quality.all)

virus.taxo %>% filter(Quality == "high quality draft" | Quality == "complete circular") %>%
  group_by(Date, order, family, TYPE, Depth, Quality) %>% tally() %>%
  mutate(new_family = ifelse(grepl("Myoriridae|Siphoviridae|Podoviridae|unknown|ambiguous|unassigned", family),
                             family,
                             "Other viral taxa")) %>%
  mutate(new_family_order = ifelse(grepl("Other", new_family),
                                   "Other taxa",
                                   paste0(order,"(",family,")"))) %>%
  ggplot(aes(x=mdy(Date), y=n, fill=family))+
  geom_bar(color="black", stat="identity", position="fill")+
  #geom_col(color="black")+
  facet_grid(Depth ~ TYPE + Quality)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())+
  ylab("Phages count")+
  xlab("Sample date")+
  ggtitle("High quality draft OR circular genomes only")

ggsave("Figures/Phage_taxo_high_qual_or_circular_only_percentage.pdf", width=11, height = 8.5)

virus.taxo %>% filter(Quality == "high quality draft" | Quality == "complete circular") %>%
  group_by(Date, order, family, TYPE, Depth, Quality) %>% tally() %>%
  mutate(new_family = ifelse(grepl("Myoriridae|Siphoviridae|Podoviridae|unknown|ambiguous|unassigned", family),
                             family,
                             "Other viral taxa")) %>%
  mutate(new_family_order = ifelse(grepl("Other", new_family),
                                   "Other taxa",
                                   paste0(order,"(",family,")"))) %>%
  ggplot(aes(x=mdy(Date), y=n, fill=family))+
  geom_bar(color="black", stat="identity", position="fill")+
  #geom_col(color="black")+
  facet_grid(Depth ~ TYPE)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())+
  ylab("Phages count")+
  xlab("Sample date")+
  ggtitle("High quality draft OR circular genomes only")

ggsave("Figures/Phage_taxo_high_qual_or_circular_only_version2_percentage.pdf", width=11, height = 8.5)

virus.taxo %>% filter(Quality == "high quality draft" | Quality == "complete circular") %>%
  group_by(Date, order, family, TYPE, Depth, Quality) %>% tally() %>%
  mutate(new_family = ifelse(grepl("Myoriridae|Siphoviridae|Podoviridae|unknown|ambiguous|unassigned", family),
                             family,
                             "Other viral taxa")) %>%
  mutate(new_family_order = ifelse(grepl("Other", new_family),
                                   "Other taxa",
                                   paste0(order,"(",family,")"))) %>%
  ggplot(aes(x=mdy(Date), y=n, fill=family))+
  geom_col(color="black")+
  facet_grid(Depth + TYPE ~. )+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())+
  ylab("Phages count")+
  xlab("Sample date")+
  ggtitle("High quality draft OR circular genomes only")

# Family and subfamily:
virus.taxo %>% filter(Quality == "high quality draft" | Quality == "complete circular") %>%
  mutate(Family_subfamily = paste(family, subfamily, sep=" ")) %>%
  group_by(Date, order, Family_subfamily, TYPE, Depth, Quality) %>% tally()  %>%
  ggplot(aes(x=mdy(Date), y=n, fill=Family_subfamily))+
  geom_bar(color="black", stat="identity", position="fill")+
  #geom_col(col="black")+
  facet_grid(Depth ~ TYPE)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())+
  ylab("Phages count")+
  xlab("Sample date")+
  ggtitle("High quality draft OR circular genomes only")+
  scale_fill_manual(values=c("grey",
                             "#FFD033","#E4B721",
                             "#BB4F1B",
                             "#05A846",
                             "#2477B0","#699DC1","#6FBBEE", #Myo
                             "#B023BC","#802788","#9D66A2",
                             "#D1003C","#E56086",
                             "#949494","#3A3A3A"
                             ))


# Family level according to Ben's suggestions.
virus.taxo %>% filter(Quality == "high quality draft" | Quality == "complete circular") %>%
  mutate(Family_subfamily = paste(family, subfamily, sep=" ")) %>%
  group_by(Date, order, family, TYPE, Depth, Quality) %>% tally()  %>%
  ggplot(aes(x=mdy(Date), y=n, fill=family))+
  geom_bar(color="black", stat="identity", position="fill")+
  #geom_col(col="black")+
  facet_grid(Depth ~ TYPE)+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank())+
  ylab("Phages count")+
  xlab("Sample date")+
  ggtitle("High quality draft OR circular genomes only")


saveRDS(virus.taxo, file = "Data/virusTaxo.rds")


ggsave("Figures/Phage_taxo_family_subfamily_colored.pdf", width = 11, height = 7)

