# Patricia Tran
# Basic plots for MAG and Phage abundance and expression
# Includes script to make the ranked plots (Figure 2)

# load packages:
library(tidyverse)
library(lubridate)

# load datasets:
# Bacteria and archaea
mag.abund <- read.table("Data/MAG_abundance_matrix_rel_abund.tsv", sep="\t", header=TRUE,check.names = FALSE) %>%
  gather("Sample.Name","relative_abund",2:17)

mag.expression <- read.table("Data/MAG_RPKM_normalized.tsv", sep="\t",
                             header=TRUE, check.names = FALSE) %>%
  gather("Sample.Name","rpkm_norm", 4:19)

#Phages
phage.abund <- read.table("Data/phages_abundance_matrix_metagenomes.tsv", sep="\t", header=TRUE,check.names = FALSE) %>%
  gather("Sample.Name","relative_abund",2:17)
phage.expression <- read.table("Data/phages_rpkm_normalized_matrix.tsv", sep="\t", header=TRUE,check.names = FALSE) %>%
  gather("Sample.Name","rpkm_norm",2:17)

####  Add metadata ####
taxonomy <- read.table("Data/MAG_taxonomy.tsv", sep="\t", header=TRUE)
sample.metadata <- read.table("Data/Samples-mendota.csv", sep=",", header=TRUE)
phage.taxonomy <- read.table("Data/Cat_phages_all_metaG_viromes.VIVID.virus-taxonomy.tsv", sep="\t", header=TRUE)
colnames(phage.taxonomy)[1] <- "Phage"

# Add metadata to MAGs
mag.abund <- mag.abund %>% left_join(taxonomy) %>% left_join(sample.metadata)
mag.expression <- mag.expression %>% left_join(taxonomy) %>% left_join(sample.metadata)

# Add metadata to phages:
phage.abund <- phage.abund %>% left_join(phage.taxonomy) %>% left_join(sample.metadata)
phage.expression <- phage.expression %>% left_join(phage.taxonomy) %>% left_join(sample.metadata)

#### PLOTS COLORS AND LAYOUT $####
cols2 <- c("p__Verrucomicrobiota"="#AF1F80",
           "p__Sumerlaeota"="grey",
           "p__Planctomycetota"="#FF0000",
           "p__Patescibacteria"="grey",
           "p__Nanoarchaeota"="grey",
           "p__Myxococcota"="grey",
           "p__Krumholzibacteriota"="grey",
           "p__Hydrogenedentota"="grey",
           "p__Gemmatimonadota"="grey",
           "p__Firmicutes"="#2a3e87",
           "p__Firmicutes_A"="#2a3e87",
           "p__Firestonebacteria"="grey",
           "p__Fibrobacterota"="grey",
           "p__Desulfobacterota"="#8242b3",
           "p__Desulfobacterota_F"="#8242b3",
           "p__Dependentiae"="grey",
           "p__Cyanobacteria"="#8BAF7B",
           "p__Chloroflexota"="#6100FF",
           "p__Chlamydiota"="grey",
           "p__Campylobacterota"="#e0540d",
           "p__Bdellovibrionota"="grey",
           "p__Bdellovibrionota_C"="grey",
           "p__Bacteroidota"="grey",
           "p__Armatimonadota"="grey",
           "p__Actinobacteriota"="#DEA955",
           "p__Acidobacteriota"="grey",
           "o__Thiomicrospirales"="#4d0c5e",
           "o__Steroidobacterales"="grey",
           "o__Pseudomonadales"="grey",
           "o__Methylococcales"="#2F83D4",
           "o__GCA-2729495"="#2FB5D4",
           "o__Burkholderiales"="#2BCFF4",
           "c__Alphaproteobacteria"="#EA9999")

#### RENAME TAXA in ABUNDANCE FILE ####
mag.abund.no.unmapped <- mag.abund %>% filter(Phylum != "unmapped")
head(mag.abund.no.unmapped)
mag.abund.no.unmapped <- mag.abund.no.unmapped %>% mutate(Taxa_to_plot = ifelse(Phylum == "p__Proteobacteria",
                                                                      yes=Class,
                                                                      no=Phylum))


mag.abund.no.unmapped <- mag.abund.no.unmapped %>% mutate(Taxa_to_plot2 = ifelse(Taxa_to_plot == "c__Gammaproteobacteria",yes=Order,
                                                                       no=Taxa_to_plot))

## Basic plots:
sample.metadata$Oxygen.level <- factor(sample.metadata$Oxygen.level, levels=c("oxic","oxycline","anoxic"))
sample.metadata$Mixing.regime <- factor(sample.metadata$Mixing.regime, 
                                        levels=c("stratified","mixed"))
ordered <- sample.metadata %>% select(Mixing.regime, 
                                      Oxygen.level,
                                      Date.YYYY.MM.DD,
                                      Sample.Depth.meters, Sample.Name) %>%
  arrange(desc(Mixing.regime), Oxygen.level) %>% distinct()

ordered
unique(ordered$Sample.Name)

mag.abund$Sample.Name <- factor(mag.abund$Sample.Name, levels=unique(ordered$Sample.Name))


mag.abund %>% 
  ggplot(aes(x=Sample.Name, y=relative_abund, fill=Phylum))+
  geom_col()+
  theme_bw()+
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle=90),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        legend.position = "bottom")+
  facet_grid(.~Mixing.regime+Oxygen.level, scale="free")+
  scale_fill_manual(values= cols2)

ggsave("Figures/Mag.abund.figure.png", width = 11, height = 8.5)



test1 <- mag.abund.no.unmapped %>% group_by(Date.YYYY.MM.DD, 
                                            Sample.Depth.meters, 
                                            Taxa_to_plot2) %>% 
  mutate(sumRelAbun = sum(relative_abund))

# Most abundant taxa:
test1 %>% select(Taxa_to_plot2, sumRelAbun) %>% distinct() %>% arrange(desc(sumRelAbun))
unique(test1$Taxa_to_plot2)


plot1.rel.abund <- test1 %>% select(Date.YYYY.MM.DD, Sample.Depth.meters, sumRelAbun, Taxa_to_plot2) %>% distinct()%>%
  ggplot(aes(x=Date.YYYY.MM.DD, y=sumRelAbun, fill=Taxa_to_plot2))+
  geom_col(position="fill")+
  scale_fill_manual(values=cols2)+
  scale_y_continuous(name = "Relative abundance (%)", 
                     labels = scales::label_percent())+
  theme_bw()+
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle=90),
        strip.background = element_blank(),
        strip.text.y = element_text(angle=0),
        panel.grid = element_blank(),
        text=element_text(color="black"))+
  facet_grid(Sample.Depth.meters~., scale="free")+
  xlab("Date")

plot1.rel.abund

ggsave(plot=plot1.rel.abund, 
       filename ="Figures/Abund.plot.pdf", width = 11, height = 8.5)

#library(plotly)
#ggplotly(plot1.rel.abund)

####  EXPRESSION ####
#### RENAME TAXA IN EXPRESSION TABLE ####
mag.expression.no.unmapped <- mag.expression %>% filter(Phylum != "unmapped")

mag.expression.no.unmapped <- mag.expression.no.unmapped %>% mutate(Taxa_to_plot = ifelse(Phylum == "p__Proteobacteria",
                                                                                yes=Class,
                                                                                no=Phylum))


mag.expression.no.unmapped <- mag.expression.no.unmapped %>% mutate(Taxa_to_plot2 = ifelse(Taxa_to_plot == "c__Gammaproteobacteria",yes=Order,
                                                                                 no=Taxa_to_plot))

head(mag.expression.no.unmapped)

#### PLOTTING EXPRESSION DATA ####

test2 <- mag.expression.no.unmapped %>% group_by(Date.YYYY.MM.DD, Sample.Depth.meters, Taxa_to_plot2) %>% 
  mutate(sum_rpkm_norm = sum(rpkm_norm))




plot2.rpkm <- test2 %>% 
  ggplot(aes(x=Date.YYYY.MM.DD, 
             y=sum_rpkm_norm, 
             fill=Taxa_to_plot2))+
  geom_col(position="fill")+
  theme_bw()+
  scale_fill_manual(values=cols2)+
  theme(panel.background = element_blank(),
        axis.text.x = element_text(angle=90),
        strip.text.y = element_text(angle=0),
        strip.background = element_blank(),
        panel.grid = element_blank())+
  facet_grid(Sample.Depth.meters~., scale="free")+
  scale_y_continuous(name = "Expression RPKM (normalized)")+
  xlab("Date")

plot2.rpkm

ggsave(plot=plot2.rpkm, 
       filename ="Figures/Expression.plot.pdf", width = 11, height = 8.5)

library(ggpubr)

ggarrange(plot1.rel.abund, plot2.rpkm, labels=c("A","B"), 
          common.legend=TRUE)

ggsave("Figures/Panels.MAGS.abund.expression.pdf", 
       width = 11, 
       height = 8.5)



#### DIFFERENT TYPES OF PLOTS TO RELATE RPKM TO EXPRESSION ####

## RPKM vs Expression:
mag.abund.expression.join <- full_join(test1, test2)

ggplot(mag.abund.expression.join, aes(x=sumRelAbun, 
                                      y=log(sum_rpkm_norm), 
                                      col=Taxa_to_plot2))+
  geom_point(size=3, alpha=0.9)+
  facet_grid(Sample.Depth.meters~Date.YYYY.MM.DD, scales = "free")+
  scale_color_manual(values=cols2)+
  scale_y_continuous(name = "log(Expression RPKM (normalized))")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y =element_text(angle=0),
        panel.grid = element_blank(),
        legend.position = "bottom")+
  xlab("Relative abundance of taxa (%)")

ggsave("Figures/RPKM_norm_vs_abund.pdf", width = 10, height = 8.5)

#### RANKING PLOTS ####

library(lubridate)
test3<- mag.abund.expression.join %>% mutate(RPKMAbundRatio = rpkm_norm/relative_abund) 

test4 <- test3 %>% select(Genome, Phylum:Species, RPKMAbundRatio, Oxygen.level, rpkm_norm, relative_abund) %>% distinct() %>% 
  mutate(sampledepth = paste0(Date.YYYY.MM.DD, "_",Sample.Depth.meters))

test4$Oxygen.level <- ifelse(test4$Date.YYYY.MM.DD == ymd("2020-10-19"),
                             "oxic mixed",test4$Oxygen.level)

unique(test4$sampledepth)
unique(test4$Oxygen.level)

ranked.taxa <- data.frame()

for (i in 1:length(unique(test4$sampledepth))){
  my_subset <- test4 %>% filter(sampledepth == unique(test4$sampledepth)[i]) %>% 
    arrange(desc(RPKMAbundRatio))
  
  ## Rank the RPKM:Abundance ratio
  rank = 1
  index = 1
  my_subset$rank <- NA
  while (my_subset$RPKMAbundRatio[index] > 0){
    my_subset$rank[index] <- rank
    rank = rank+1
    index = index +1
  }
  
  my_subset$rank_cat <- NA
  my_subset$rank_cat[(index-1-10):(index-1)] <- "Bottom 10"
  my_subset$rank_cat[1:10] <- "Top 10"
  
  my_subset$rank_cat_to_plot <- NA
  my_subset$rank_cat_to_plot[(index-10):(index-1)] <- seq(11,20,1)
  my_subset$rank_cat_to_plot[1:10] <- seq(1,10,1) 
  
  ## Rank by Relative abundance:
  my_subset <- my_subset %>% arrange(desc(relative_abund))
  
  rank.abund = 1
  index.abund = 1
  my_subset$rank_relabund <- NA
  while (my_subset$relative_abund[index.abund] > 0){
    my_subset$rank_relabund[index.abund] <- rank.abund
    rank.abund = rank.abund+1
    index.abund = index.abund +1
  }
  
  rank.abund
  index.abund

  
  my_subset$rank_cat_abund <- NA
  my_subset$rank_cat_abund[(index.abund-10):(index.abund-1)] <- "Bottom 10"
  my_subset$rank_cat_abund[1:10] <- "Top 10"
  
  my_subset$rank_cat_to_plot_abund <- NA
  my_subset$rank_cat_to_plot_abund[(index.abund-10):(index.abund-1)] <- seq(11,20,1)
  my_subset$rank_cat_to_plot_abund[1:10] <- seq(1,10,1) 
  
  
  # Rank by expression
  my_subset <- my_subset %>% arrange(desc(rpkm_norm))
  
  rank.expression = 1
  index.expression = 1
  my_subset$rank_expression <- NA
  while (my_subset$rpkm_norm[index.expression] > 0){
    my_subset$rank_expression[index.expression] <- rank.expression
    rank.expression = rank.expression+1
    index.expression = index.expression +1
  }
  
  my_subset$rank_cat_expression <- NA
  my_subset$rank_cat_expression[(index.expression-10):(index.expression-1)] <- "Bottom 10"
  my_subset$rank_cat_expression[1:10] <- "Top 10"
  
  
  my_subset$rank_cat_to_plot_expression <- NA
  my_subset$rank_cat_to_plot_expression[(index.expression-10):(index.expression-1)] <- seq(11,20,1)
  my_subset$rank_cat_to_plot_expression[1:10] <- seq(1,10,1) 
  
  
  ranked.taxa <- rbind(ranked.taxa, my_subset)
}



ordered <- ranked.taxa %>% select(Date.YYYY.MM.DD, Sample.Depth.meters, sampledepth) %>%
  arrange(desc(Date.YYYY.MM.DD), desc(Sample.Depth.meters)) %>% distinct()

ordered
unique(ordered$sampledepth)

# Reorder factors:
ranked.taxa$rank_cat <- factor(ranked.taxa$rank_cat, levels=c("Top 10","Bottom 10"))
ranked.taxa$rank_cat_expression <- factor(ranked.taxa$rank_cat_expression, levels=c("Top 10","Bottom 10"))
ranked.taxa$rank_cat_abund <- factor(ranked.taxa$rank_cat_abund, levels=c("Top 10","Bottom 10"))

ranked.taxa$Oxygen.level <- factor(ranked.taxa$Oxygen.level, levels=c("oxic","oxycline","anoxic","oxic mixed"))
ranked.taxa$sampledepth <- factor(ranked.taxa$sampledepth, 
                                 levels = unique(ordered$sampledepth))


my_rank_theme <-   theme_bw()+
  theme(axis.text.x = element_text(angle=0, hjust=1),
        panel.grid  = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(angle=0),
        strip.text.y = element_text(angle=0),
        legend.position="bottom",
        panel.border = element_blank(),
        axis.ticks = element_blank())

# Make a ranked plot for Ratio, abundance and Expression
rank.ratio.plot <- ggplot(ranked.taxa, aes(x=rank, y=sampledepth, fill=Taxa_to_plot2, color=rank_cat))+
  geom_tile(size=0.5)+
  scale_fill_manual(values=cols2)+
  scale_color_manual(values=c("black","black"),
                     na.value = NA)+
  my_rank_theme+
  facet_grid(Oxygen.level~., scales = "free", space = "free")+
  ggtitle("Ranked by ratio of activity:abundance")

rank.ratio.plot 

ggsave(plot= rank.ratio.plot,
       filename = "Figures/Top_Bottom_RPKM_abund_ratio.pdf", width = 11, height = 6)

topbottom_rank.plot <- ggplot(ranked.taxa %>% filter(!is.na(rank_cat)), 
       aes(x=rank_cat_to_plot, y=sampledepth, fill=Taxa_to_plot2))+
  geom_tile(size=1, color="white")+
  scale_fill_manual(values=cols2)+
  scale_color_manual(values=c("black","black"),
                     na.value = NA)+
  my_rank_theme+
  facet_grid(Oxygen.level~rank_cat, scales = "free", space = "free")+
  ggtitle("Top 10 and Bottom 10 ranked ratios")

topbottom_rank.plot


# Ranked Abundance:
unique(ranked.taxa$rank_cat_abund)
unique(ranked.taxa$rank_relabund)

rank.abund.plot <- ggplot(ranked.taxa, aes(x=rank_relabund, 
                                           y=sampledepth, 
                                           fill=Taxa_to_plot2, color=rank_cat_abund))+
  geom_tile(size=0.5)+
  scale_fill_manual(values=cols2)+
  scale_color_manual(values=c("black","black"),
                     na.value = NA)+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        panel.grid  = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(angle=0),
        strip.text.y = element_text(angle=0),
        legend.position="bottom",
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  facet_grid(Oxygen.level~., scales = "free", space = "free")+
  ggtitle("Ranked by relative abundance")

rank.abund.plot

ggsave(plot = rank.abund.plot,
       filename = "Figures/ranked_abund.plot.pdf", width = 11, height = 6)

topbottom_rank.abund.plot <- ggplot(ranked.taxa %>% filter(!is.na(rank_cat_abund)), 
                              aes(x=rank_cat_to_plot_abund, 
                                  y=sampledepth, fill=Taxa_to_plot2))+
  geom_tile(size=1, color="white")+
  scale_fill_manual(values=cols2)+
  scale_color_manual(values=c("black","black"),
                     na.value = NA)+
  my_rank_theme+
  facet_grid(Oxygen.level~rank_cat_abund, scales = "free", space = "free")+
  ggtitle("Top 10 and Bottom 10 abundance")

topbottom_rank.abund.plot

# Rank expression:
unique(ranked.taxa$rank_expression)

rank.expression.plot <- ggplot(ranked.taxa, aes(x=rank_expression, y=sampledepth, fill=Taxa_to_plot2, color=rank_cat_expression))+
  geom_tile(size=0.5)+
  scale_fill_manual(values=cols2)+
  scale_color_manual(values=c("black","black"),
                     na.value = "white")+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1),
        panel.grid  = element_blank(),
        strip.background = element_blank(),
        strip.text.x = element_text(angle=0),
        strip.text.y = element_text(angle=0),
        legend.position="bottom",
        panel.border = element_blank(),
        axis.ticks = element_blank())+
  facet_grid(Oxygen.level~., scales = "free", space = "free")+
  ggtitle("Ranked by expression")

rank.expression.plot

ggsave(plot = rank.expression.plot, filename= "Figures/ranked_expression.pdf",
       width = 11, height = 6)


topbottom_rank.expression.plot <- ggplot(ranked.taxa %>% filter(!is.na(rank_cat_expression)), 
                                    aes(x=rank_cat_to_plot_expression, 
                                        y=sampledepth, fill=Taxa_to_plot2))+
  geom_tile(size=1, color="white")+
  scale_fill_manual(values=cols2)+
  scale_color_manual(values=c("black","black"),
                     na.value = NA)+
  my_rank_theme+
  facet_grid(Oxygen.level~rank_cat_expression, scales = "free", space = "free")+
  ggtitle("Top 10 and Bottom 10 expression")

topbottom_rank.expression.plot


ggarrange(rank.abund.plot, rank.expression.plot, rank.ratio.plot, nrow=3,
          common.legend=TRUE, labels=c("A","B", "C"))

ggsave("Figures/Ranked_all.pdf", width = 11, height = 8.5)

ggarrange(topbottom_rank.abund.plot, topbottom_rank.expression.plot, topbottom_rank.plot, nrow=3,
          common.legend=TRUE, labels=c("A","B", "C"))

ggsave("Figures/ranked_top10bottom10.pdf", width = 8.5, height = 11)

### Panel with top only:
ranked.taxa.top <- ranked.taxa %>% 
  filter(!is.na(rank_cat_expression) | !is.na(rank_cat_abund) | !is.na(rank_cat))

phages.host <- read.csv("Data/Host_prediction_to_genome_m90.csv")

ranked.taxa.top <- ranked.taxa.top %>% mutate(hasPhage = ifelse(Genome %in% phages.host$Host.genome,
                                             yes="*",
                                             no=NA))

#View(ranked.taxa.top)

ranked.taxa.top.1 <- ranked.taxa.top %>% select(Taxa_to_plot2, sampledepth, Oxygen.level,
                           rank_cat_expression, rank_cat_to_plot_expression,
                           rank_cat_abund, rank_cat_to_plot_abund,
                           rank_cat, rank_cat_to_plot,
                           hasPhage) %>%
  gather("rank","rank_to_plot", c(rank_cat_to_plot_expression,
                                  rank_cat_to_plot_abund,
                                  rank_cat_to_plot)) %>%
  filter(rank_to_plot <= 10)



ranked.taxa.top.1$rank <- factor(ranked.taxa.top.1$rank,
                                  levels=c("rank_cat_to_plot_abund",
                                           "rank_cat_to_plot_expression",
                                           "rank_cat_to_plot"))


ggplot(ranked.taxa.top.1,
  aes(x=rank_to_plot, 
      y=sampledepth, 
      fill=Taxa_to_plot2))+
  geom_tile(size=1, color="white")+
  geom_text(aes(label=hasPhage), 
            size =6)+
  scale_fill_manual(values=cols2)+
  #scale_color_manual(values="black",
  #                   na.value = NA)+
  my_rank_theme+
  facet_grid(Oxygen.level~rank, scales = "free", space = "free")+
  ggtitle("Top 10")

# This is Figure 2B:
ggsave("Figures/Top10-paneled.rank.pdf", width = 11, height = 6)


#### PHAGES - are there  any phages expressed at every time point?####

phage.exp.pattern <- phage.expression %>% filter(rpkm_norm > 0) %>%
  group_by(Phage) %>% tally() %>% arrange(desc(n)) %>%
  mutate(Dataset = "Expressed in")

phage.abund.pattern <- phage.abund %>% filter(relative_abund > 0) %>%
  group_by(Phage) %>% tally() %>% arrange(desc(n)) %>%
  mutate(Dataset = "Found in")

## No, at most a phage was expressed in 9 samples.
## Therer are mutiple phage abundant in multiple samples

phage.abund.expression.pattern <- full_join(phage.exp.pattern, phage.abund.pattern) %>%
  spread(Dataset, n) %>%
  select(Phage, `Found in`, `Expressed in`) %>%
  arrange(desc(`Found in`), desc(`Expressed in`))


head(phages.host)

phages.abund.expression.with.host.information <- full_join(phage.abund.expression.pattern, 
                                                           phages.host, 
                                                           by=c("Phage"="Virus"))

phages.abund.expression.with.host.information <- phages.abund.expression.with.host.information %>% mutate("CyanoHost" = ifelse(grepl("Cyano", 
                                                                                    Host.taxonomy),
                                                                              yes= "Cyanobacteria host",
                                                                              no = "Other bacterial host"))

#View(phages.abund.expression.with.host.information)

phages.host



#View(phage.abund.expression.pattern)

phage.abund.expression.pattern.tally <- phage.abund.expression.pattern %>%
  group_by(`Found in`,`Expressed in`) %>% tally()

phage.abund.expression.pattern.tally$`Expressed in` <- replace_na(phage.abund.expression.pattern.tally$`Expressed in`,
                                                                  0)

phage.abund.expression.pattern.tally <- phages.abund.expression.with.host.information %>% 
  select(Phage, `Found in`,`Expressed in`, CyanoHost) %>% distinct() %>%
  group_by(`Found in`,`Expressed in`, CyanoHost)%>%
  tally()

phage.abund.expression.pattern.tally$`Expressed in` <- replace_na(phage.abund.expression.pattern.tally$`Expressed in`,
                                                                  0)
  
head(phage.abund.expression.pattern.tally)
phage.abund.expression.pattern.tally %>%
  ggplot(aes(x=`Found in`, 
             y=`Expressed in`, 
             fill=n))+
  geom_tile(col="black", size=1)+
  scale_x_continuous(breaks= c(0:16))+
  scale_y_continuous(breaks= c(0:16),
                     limits = c(0,16))+
  theme_bw()+
  facet_wrap(CyanoHost~.) +
  theme(panel.grid = element_blank(),
        strip.background = element_blank())+
  geom_text(aes(label=n), col="dark grey", size=4)+
  xlab("Phages found in how many samples  (max = 16)")+
  ylab("Phages active in how many samples (max = 16)")
  #scale_color_manual(values=c("dark green","black"))

# Sup Figure 5
ggsave("Figures/PhagesFoundInAndExpressedInHowManySamples.pdf",
       width = 10,
       height = 5)  
  
## Add taxonomy information
phages.abund.expression.with.host.information <- phages.abund.expression.with.host.information %>% separate(Host.taxonomy,
                                                           into=c("domain","phylum","class",
                                                                  "order","genus","family","species"),
                                                           sep=";",
                                                           remove=FALSE)

taxa_phages_patterns <- phages.abund.expression.with.host.information %>% select(`Found in`,`Expressed in`,
                                                         class, order, genus, family) %>%  distinct()

taxa_phages_patterns$`Expressed in` <- replace_na(taxa_phages_patterns$`Expressed in`,0)

count_info <- phage.abund.expression.pattern %>% 
  group_by(`Found in`,`Expressed in`) %>% tally()

count_info$`Expressed in` <- replace_na(count_info$`Expressed in`, 0)



taxa_phages_patterns_count <- left_join(taxa_phages_patterns, count_info)
colnames(taxa_phages_patterns_count)[7] <- "PhagesCountMatching"

#View(taxa_phages_patterns_count)

ggplot(taxa_phages_patterns_count, aes(x=`Found in`,
                                       y=`Expressed in`,
                                       size=PhagesCountMatching,
                                       col=class))+
  geom_point()

# get top expressed phage
phages.with.MAGhost.expressed <- phages.abund.expression.with.host.information %>% filter(grepl("Ga",Host.genome)) %>%
  filter(!is.na(`Expressed in`))
head(phage.expression)

phage.expression.subset <- phage.expression %>% filter(Phage %in% phages.with.MAGhost.expressed$Phage)

# How many times is this host  expressed?

head(phage.expression)
head(mag.abund)

list.unique.phages.with.mag.host <- unique(phage.expression.subset$Phage)

phages.with.MAGhost.expressed <- phages.with.MAGhost.expressed %>% 
  left_join(mag.expression, by=c("Host.genome"="Genome"))

phages.with.MAGhost.expressed2 <- phages.with.MAGhost.expressed %>% 
  filter(rpkm_norm > 0) %>%
  select(Phage, `Found in`, `Expressed in`, rpkm_norm) %>%
  group_by(Phage, `Found in`, `Expressed in`) %>% tally() %>%
  arrange(desc(`Found in`),desc(`Expressed in`), desc(n)) %>%
  mutate(HostExpressedInThisManySamples = n) %>%
  select(-n) %>% 
  left_join(phages.abund.expression.with.host.information) %>%
  filter(grepl("Ga", Host.genome))

write.table(phages.with.MAGhost.expressed2, "Data/Phages_with_Host_both_expressed.tsv",
            row.names = FALSE,
            sep="\t",
            quote = FALSE)


## Expression patterns of the 55 MAGS:
mags55.ranked <- ranked.taxa %>% filter(Genome %in% phages.host$Host.genome)
mags55.ranked

# how many of those 55  are not expressed.
mags55.ranked %>% group_by(Genome) %>% mutate(SumRPKM = sum(rpkm_norm)) %>%
  select(Genome, SumRPKM)  %>% distinct() %>% filter(SumRPKM == 0) %>% nrow()

test1 <- mag.abund.expression.join %>% filter(Genome %in% phages.host$Host.genome)

colnames(test1)
max.value <- max(test1$rpkm_norm)

test1 <- test1[,c(1,2,26)]

test1 <- test1 %>% spread(Sample.Name, rpkm_norm)

colnames(test1)
data <- as.matrix(test1[,2:17])
row.names(data) <- test1$Genome

write.table(data, "Data/Expression_55MAGS.tsv", row.names = TRUE,
            quote = FALSE, sep="\t")





