# CoverM on metagenomes and metatranscriptomes:
library(tidyverse)
library(lubridate)
detach("package:MASS", unload = TRUE)


coverM_metaG <- read.csv("Data/431MAGS-vs-Metagenomes-CoverM.tsv", sep="\t")
coverM_metaG <- coverM_metaG %>% gather("Sample","RelativeAbundance",2:17)

coverM_metaT <- read.csv("Data/431MAGS-vs-Metatranscriptomes-CoverM.tsv", sep="\t")
coverM_metaT <- coverM_metaT %>% gather("Sample","RelativeAbundance",2:17)

head(coverM_metaG)
head(coverM_metaT)

coverM_metaG_metaT <- full_join(coverM_metaG, coverM_metaT)

unique(coverM_metaG_metaT$Sample)
coverM_metaG_metaT$Sample <- str_replace(coverM_metaG_metaT$Sample,".fastq.Relative.*","")
unique(coverM_metaG_metaT$Sample)
coverM_metaG_metaT$Sample <- str_replace(coverM_metaG_metaT$Sample,"X5","5")
unique(coverM_metaG_metaT$Sample)


coverM_match <- read.csv("Data/coverM_match_samples.csv")
coverM_metaG_metaT <- left_join(coverM_metaG_metaT, coverM_match)

coverM_metaG_metaT$Date <- mdy(coverM_metaG_metaT$Date)

#hcgA_genes 

hgca <- read.table("Data/MAG_Characteristics_Binned_hgcA_LakeMendotaPairedViromeMetagenome-PQT2020.tsv", sep="\t", header=T)

colnames(hgca)

hgca.desulfo <- hgca %>% filter(Phylum=="p__Desulfobacterota") %>% 
  dplyr::select(V1, V2, completeness, Domain, Phylum, Class, Order, Family, Genus, Species)

plot1 <- coverM_metaG_metaT %>% filter(Genome %in% hgca.desulfo$V2) %>%
  filter(RelativeAbundance != 0) %>%
  ggplot(aes(x=Date, y=RelativeAbundance, col=Genome, shape=Type, group=interaction(Genome, Type)))+
  geom_point(size=2, alpha=0.8)+
  facet_grid(Depth~Type, scales="free")+
  theme_bw()+
  geom_line(orientation = "x", alpha=0.5, size=1.5)+
  ggtitle("Desulfobacteria MAGs which contains hgcA genes", 
          subtitle="Relative abundance calculated with coverM")+
  theme(panel.grid = element_blank(),
        text=element_text(size=12, color="black"),
        strip.background = element_blank())+
  ylab("Relative abundance (%)")+
  xlab("Date (2020)")


plot1

ggsave(plot1, filename = "Figures/Desulfobacteria_with_hcgA_remove_0.pdf", width = 8, height = 9)


# HgcA total:
coverM_all_hcgA_containing_bins <- left_join(coverM_metaG_metaT, hgca, by=c("Genome"="V2")) %>%
  filter(!is.na(Marker.lineage))

plot2 <- coverM_all_hcgA_containing_bins %>%
  filter(RelativeAbundance != 0) %>%
  ggplot(aes(x=Date, y=RelativeAbundance, col=Phylum, shape=Type, group=interaction(Genome, Type)))+
  geom_point(size=3, alpha=0.8)+
  facet_grid(Depth~Type, scales="free")+
  theme_bw()+
  geom_line(orientation = "x", alpha=0.5, size=1.5)+
  ggtitle("All MAGs which contains hgcA genes", 
          subtitle="Relative abundance calculated with coverM")+
  theme(panel.grid = element_blank(),
        text=element_text(size=12, color="black"),
        strip.background = element_blank())+
  ylab("Relative abundance (%)")+
  xlab("Date (2020)")

plot2
ggsave(plot2 , filename = "Figures/AllMAGs_with_hcgA_remove_0.pdf", width = 8, height = 9)

library(ggpubr)
ggarrange(plot2,plot1, labels=c("A","B"))

ggsave("Figures/Panels_hgcA_remove0.pdf", width = 13, height = 8)

## CoverM Bar Plot

mag_char <- read.table("Data/MAG_Characteristics.tsv", sep="\t", header=T)

mag_char <- mag_char %>% select(genome, Domain, Phylum)
coverM_metaG_metaT <- left_join(coverM_metaG_metaT, mag_char, by=c("Genome"="genome"))

coverM_metaG_metaT$Phylum <- str_replace(coverM_metaG_metaT$Phylum, "p__","")


ggplot(coverM_metaG_metaT %>% filter(!is.na(Phylum)), aes(x=SampleName, y=RelativeAbundance, fill=Phylum)) +
  geom_col()+
  facet_grid(.~Type, scales="free" )+
  theme_bw()+
  scale_fill_manual(breaks = c("Bdellovibrionota","Bdellovibrionota_C","Campylobacterota","Chlamydiota","Dependentiae","Fibrobacterota","Firestonebacteria","Gemmatimonadota","Hydrogenedentota","Krumholzibacteriota","Nanoarchaeota","Sumerlaeota","Armatimonadota","Chloroflexota","Firmicutes","Firmicutes_A","Myxococcota","Patescibacteria","Planctomycetota",
                               "Acidobacteriota","Verrucomicrobiota","Desulfobacterota","Desulfobacterota_F","Bacteroidota","Gammaproteobacteria","Cyanobacteria","Alphaproteobacteria","Actinobacteriota"),
                    values=c("#dbdbdb","#dbdbdb","#dbdbdb","#dbdbdb","#dbdbdb","#dbdbdb","#dbdbdb","#dbdbdb","#dbdbdb","#dbdbdb","#dbdbdb","#dbdbdb","#dbdbdb","#dbdbdb","#dbdbdb","#dbdbdb","#dbdbdb","#dbdbdb","#dbdbdb",
                             "#cc903b","#c22151","#a966d9","#a966d9","#693705","#5e95bf","#099e3e","#f748b1","#e8b061"))+
  theme(panel.grid = element_blank(),
        text=element_text(size=12, color="black"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle=90))

library(plotly)
#ggplotly()

ggsave("Figures/coverM_bar_plots_no_unbinned.pdf", width = 11, height = 8)


### COVERM on Phages:
# Run MMSEQ on phages from vibrant, then run coverm on the mmseq phages with the viromes, metagenomes and metatranscriptomes:


phages_MG_coverM <- read.table("Data/MMSEQ_phages-vs-Metagenomes-CoverM.tsv", sep="\t", header=TRUE) %>% gather("Sample","RelativeAbundance",2:17)
phages_MT_coverM <- read.table("Data/MMSEQ_phages-vs-Metatranscriptomes-CoverM.tsv", sep="\t", header=TRUE)%>% gather("Sample","RelativeAbundance",2:17)
phages_Viromes_coverM <- read.table("Data/MMSEQ_phages-vs-Viromes-CoverM.tsv", sep="\t", header=TRUE)%>% gather("Sample","RelativeAbundance",2:15)

coverM_phages <- full_join(phages_MG_coverM, phages_MT_coverM)
coverM_phages <- full_join(coverM_phages, phages_Viromes_coverM)

coverM_phages$Sample <- str_replace(coverM_phages$Sample,".fastq.Relative.*","")
coverM_phages$Sample <- str_replace(coverM_phages$Sample,"X5","5")
coverM_phages$Sample <- str_replace(coverM_phages$Sample,".fastq.*","")
unique(coverM_phages$Sample)

coverM_phages <- left_join(coverM_phages, coverM_match)

ggplot(coverM_phages, aes(x=Type, y=RelativeAbundance, fill=Genome))+
  geom_col(stat="identity")+
  facet_grid(Depth~mdy(Date), scales="free")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        text=element_text(size=12, color="black"),
        strip.background = element_blank(),
        axis.text.x = element_text(angle=90))+
  ylab("Relative abundance in sample (%)")+
  xlab("Sample Type")+
  scale_fill_manual(values=c("red","grey"),
                    labels = c("Dereplicated phages (all)", "Unmapped"))+
  ggtitle("VIBRANT on all samples > MMSEQ to cluster > Get unique phages > CoverM")

ggsave("Figures/CoverM_Phages_vs_mg_mt_viromes.pdf", width = 11, height = 8.5)

head(coverM_phages)

ggplot(coverM_phages %>% filter(Genome != "unmapped" & !is.na(RelativeAbundance)) , 
       aes(x=mdy(Date), y=RelativeAbundance, col=as.factor(Depth), shape=Type))+
  geom_point(size=3)+
  geom_line(orientation = "x")+
  facet_wrap(Type~.)+
  ylab("Relative Abundance in sample (%)")+
  xlab("Date")+
  scale_x_date(labels=mdy(unique(coverM_phages$Date)),
               breaks=mdy(unique(coverM_phages$Date)))+
  theme_bw()+
  theme(panel.grid.minor= element_blank(),
        axis.text.x = element_text(angle=90),
        strip.background = element_blank())+
  scale_color_manual(values=c("yellow","orange","blue","purple"),
                     labels=c("5m","10m","15m","23.5m"))

ggsave("Figures/CoverM_phages_over_time_depth_noMT.pdf", width = 11, height = 4)
