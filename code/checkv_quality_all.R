# CheckV quality of all phages mt and mg
library(data.table)
library(tidyverse)

quality.checkv <- fread("Data/quality_summary-2.tsv", sep="\t", header=TRUE)
completeness.checkv <- fread("Data/completeness-1.tsv", sep="\t", header=TRUE)


checkv_all <- full_join(quality.checkv, completeness.checkv, by=c("contig_id", "contig_length"))

checkv_all %>% group_by(checkv_quality) %>% tally()

phage.genome.qual.plot <- checkv_all %>% mutate(phagetype = ifelse(test = contig_length/1000 < 15, 
                                         yes="miniphage", 
                                         no=ifelse(test=contig_length/1000 > 200, 
                                                   yes="mega-phage", 
                                                   no="normal phage"))) %>%
ggplot(aes(x=checkv_quality, y=contig_length/1000))+
  geom_boxplot()+
  geom_point(position=position_jitterdodge(jitter.width=0, dodge.width = 0, seed = 1234),
             aes(color=phagetype), show.legend = F)+
  theme_bw()+
  ylab("Genome size (kbp)")+
  xlab("Viral genome quality (CheckV)")+
  geom_hline(yintercept =15000/1000)+
  geom_hline(yintercept = 20000/100)

ggsave(plot = phage.genome.qual.plot,
       filename = "Figures/check_quality_genome_size.pdf", width=8, height=8)

# Add whether it's a mini-phage, mega-phage, or other
checkv_all <- checkv_all %>% mutate(phagetype = ifelse(test = contig_length/1000 < 15, 
                                                       yes="miniphage", 
                                                       no=ifelse(test=contig_length/1000 > 200, 
                                                                 yes="mega-phage", 
                                                                 no="normal phage"))) 

checkv_all <- checkv_all %>% mutate(sample_name = str_replace(contig_id, "_.*","")) %>%
  mutate(sample_type = ifelse(test = as.numeric(str_sub(sample_name,-2))>=73,
                              yes = "Virome",
                              no = "Metagenome"))

checkv_all %>%
  group_by(sample_type, checkv_quality, phagetype) %>% tally() %>%
  ggplot() +
  geom_col(aes(x=checkv_quality, y=n, fill=phagetype),
           position="dodge")+
  theme_bw()+
  theme(panel.grid = element_blank())

checkv_all %>% filter(checkv_quality %in% c("Complete","High-quality")) %>%
  group_by(sample_type, checkv_quality, phagetype) %>% tally() %>%
  ggplot(aes(x=sample_type, y=n, fill=phagetype))+
  geom_col(aes(x=sample_type, y=n, fill=phagetype),
           position="dodge")+
  theme_bw()+
  theme(panel.grid = element_blank())+
  ylab("Number of phages (high-quality or complete only)")+
  ylab("Sample origin")
  
ggsave("Figures/good_phages_by_sample_origin.pdf", width = 5, height = 5)


# Relationship between lytic and lysogenic phages: and genome quality:
#Compare Lytic and Lysogenic Phages

lytic.MG <- fread("Data/Lytic_phages.txt", sep="\t", header=FALSE)
lyso.MG <- fread("Data/Lysogenic_phages.txt", sep="\t", header=FALSE)

lytic.virome <- fread("Data/Lytic_phages_virome.txt", sep="\t", header=FALSE)
lyso.virome <- fread("Data/Lysogenic_phages_virome.txt", sep="\t", header=FALSE)

all_lyso_lytic <- rbind(lytic.MG, lyso.MG, lytic.virome, lyso.virome)

head(all_lyso_lytic)
dim(all_lyso_lytic)    
colnames(all_lyso_lytic) <- c("contig_id","type_by_vibrant","sample_type")

dim(checkv_all)
checkv_all <- left_join(checkv_all, all_lyso_lytic)  

head(checkv_all)


pairs <- read.csv("Data/MG_virome_pairs.csv")

head(checkv_all)
head(pairs)

checkv_all <- left_join(checkv_all, pairs, by=c("sample_name" = "GaID"))  

checkv_all %>% group_by(Paired.Sample.Number, checkv_quality, type_by_vibrant, sample_type) %>% tally() %>%
  filter(Paired.Sample.Number != "Unpaired") %>% 
  ggplot(aes(x=as.numeric(Paired.Sample.Number), y=n, 
             col=sample_type, 
             group=Paired.Sample.Number))+
  geom_point(size=2)+
  geom_line(orientation="y", color="grey", alpha=0.5)+
  facet_wrap(checkv_quality~type_by_vibrant, scales="free")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        text= element_text(color="black"))+
  ylab("Count")+
  scale_color_manual(values=c("purple","orange"))+
  scale_x_continuous(breaks = seq(1,14,1))+
  xlab("Sample pair")

ggsave("Figures/Comparing_lytic_vs_lysogenic_by_type_quality_and_pair_2.pdf", width = 11, height = 6)

compare.phages.good.phages <- checkv_all %>% group_by(Paired.Sample.Number, checkv_quality, 
                        type_by_vibrant, sample_type, SampleName_no_type, Oxygen) %>% 
  tally() %>%
  mutate(Oxygen_f = factor(Oxygen, levels=c("Oxic","Oxic-Anoxic","Anoxic"))) %>%
  filter(checkv_quality %in% c("Complete","High-quality")) %>%
  ggplot(aes(y= SampleName_no_type, 
             x=n, 
             col=sample_type, 
             group=Paired.Sample.Number,
             shape= checkv_quality))+
  geom_point(size=5)+
  geom_line(color="black")+
  facet_grid(Oxygen_f~sample_type+checkv_quality, scales="free", space = "free")+
  theme_bw()+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        text= element_text(color="black"))+
  xlab("Count")+
  ylab("Sample pair")+
  #scale_y_continuous(breaks = seq(1,14,1))+
  #ggtitle("Complete & high-quality viruses only")+
  scale_color_manual(values=c("turquoise","blue"))

compare.phages.good.phages

ggsave("Figures/PairedMG_virome_complete_HQ_viruses.pdf", 
       width = 8.5, 
       height = 11)

## Add predicted host taxonomy to checkv file:

links <- read.csv("Data/Host_prediction_to_genome_m90.csv")
head(links)
head(checkv_all)

checkv_all_with_iphop <- left_join(checkv_all, links, by=c("contig_id"="Virus"))

write.table(checkv_all_with_iphop, "Data/Phages_all_qualities_and_sizes_with_taxonomy.tsv", sep="\t", row.names = FALSE,
            quote = FALSE)


checkv_all_with_iphop %>% filter(miuvig_quality == "High-quality")

write.table(checkv_all_with_iphop %>% filter(miuvig_quality == "High-quality"), 
            "Data/Phages_all_qualities_and_sizes_with_taxonomy_miuvig_highqual.tsv", sep="\t", row.names = FALSE,
            quote = FALSE)

write.table(checkv_all_with_iphop %>% filter(checkv_quality %in% c("Complete","High-quality")), 
            "Data/Phages_all_qualities_and_sizes_with_taxonomy_checkv_complete_highqual.tsv", sep="\t", row.names = FALSE,
            quote = FALSE)


checkv_to_use <- checkv_all %>% filter(checkv_quality %in% c("Complete","High-quality"))

my_theme <- theme_bw()+
  theme(panel.grid=element_blank(),
        text = element_text(size=12, color="black"))

cols_oxygen <- c("Oxic"="#48AA72",
                 "Oxic-Anoxic"="red",
                 "Anoxic"="#3E1E53")

library(ggrepel)

data.to.view <- checkv_to_use %>% group_by(sample_type, Oxygen, Paired.Sample.Number, SampleName_no_type) %>% 
  tally()
#View(data.to.view)


# Figure 3
panelA <- checkv_to_use %>% group_by(sample_type, Oxygen, Paired.Sample.Number, SampleName_no_type) %>% 
  tally()%>%
  ggplot(aes(x=sample_type, y=n, 
             group=Paired.Sample.Number, col=Oxygen))+
  geom_point(size=2)+
  geom_line()+
  my_theme+
  scale_color_manual(values = cols_oxygen)+
  geom_text_repel(aes(label=SampleName_no_type))

panelA

panelB <- checkv_to_use %>% group_by(sample_type, Oxygen, sample_name, type_by_vibrant) %>% tally() %>%
  ggplot(aes(x=sample_type, n, group=sample_type, col=Oxygen))+
    geom_violin()+
  geom_point()+
  facet_wrap(type_by_vibrant~., scale="free_y")+
  my_theme+
  theme(strip.background = element_blank())+
  stat_summary(fun = "mean", size = 1, col="red")+
  stat_summary(fun="sd", size=1, col="blue")+
  scale_color_manual(values = cols_oxygen)

panelB

plotA <- plotA + my_theme

ggarrange(plotA, panelA, panelB, labels=c("A","B","C"), nrow=1, 
          widths = c(0.4,0.4,1), common.legend = TRUE)

ggsave("Figures/Figure2_phages_ttest_types.pdf", width = 11, height = 5)


library(PairedData)

checkv_to_use %>% filter(Paired.Sample.Number != "Unpaired") %>% 
  

metagenomes<- subset(checkv_all %>% filter(Paired.Sample.Number != "Unpaired"),  
                     sample_type == "Metagenome", Putative.phages.found.percent,
                     drop = TRUE)

viromes<- subset(vibrant_results %>% filter(Paired.Sample.Number != "Unpaired"),  Type == "Virome", Putative.phages.found.percent,
                 drop = TRUE)

paired.data <- paired(metagenomes,viromes)
head(paired.data)

library(ggpubr)
ggpaired(paired.data, cond1="metagenomes", 
         cond2="viromes")+
  ylab("Percent of data that are phages")+
  xlab("Sample type")+
  my_fave_theme

ggsave("Figures/Pairs.phages.boxplot.pdf", units = "in", height = 8, width = 8)

# compute the difference
d <- with(vibrant_results %>% filter(Paired.Sample.Number != "Unpaired"), 
          Putative.phages.found.percent[Type == "Metagenome"] - Putative.phages.found.percent[Type == "Virome"])


# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = p-value = 0.8523


#From the output, the p-value is greater than the significance level 0.05 implying that the distribution of the differences (d) are not significantly different from normal distribution. In other words, we can assume the normality.
res <- t.test(metagenomes, viromes, paired = TRUE)
res







