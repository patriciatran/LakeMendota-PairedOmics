library(BiodiversityR) # also loads vegan
library(ggplot2)
library(readxl)
library(ggsci)
library(ggrepel)
library(ggforce)

BioR.theme <- theme(
  panel.background = element_blank(),
  panel.border = element_blank(),
  panel.grid = element_blank(),
  axis.line = element_line("gray25"),
  text = element_text(size = 12),
  axis.text = element_text(size = 10, colour = "gray25"),
  axis.title = element_text(size = 14, colour = "gray25"),
  legend.title = element_text(size = 14),
  legend.text = element_text(size = 14),
  legend.key = element_blank())

sp.matrix <- read.table("Data/MAG_abundance_matrix_rel_abund.tsv", 
                          sep="\t", header=TRUE)

summary(sp.matrix)

library(tidyverse)


row.names(sp.matrix) <- str_replace(sp.matrix$Genome,"X","")
sp.matrix <- sp.matrix[,-1]
sp.matrix <- t(sp.matrix)



sp.env <- read.csv("Data/Samples-mendota.csv")
row.names(sp.env) <- sp.env[,1]
sp.env <- sp.env[,-1]

head(sp.env)

#sp.env %>% filter(Date.YYYY.MM.DD == "2020-10-19")

sp.env <- sp.env %>% mutate(Mixing_Oxygen = ifelse(Date.YYYY.MM.DD == "2020-10-19", 
                                                   yes="Mixed",
                                         no = paste("Stratified", Oxygen.level)))

unique(sp.env$Mixing_Oxygen)

Ordination.model1 <- metaMDS(sp.matrix, distance='bray', 
                             k=2, trymax=1, 
                             autotransform=TRUE,
                             noshare=0.1, expand=TRUE, trace=1, plot=FALSE)

plot1 <- ordiplot(Ordination.model1, choices=c(1,2))

sites.long1 <- sites.long(plot1, env.data=sp.env)

head(sites.long1)

plotgg1 <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab("NMDS1") +
  ylab("NMDS2") +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
  geom_point(data = sites.long1, 
             aes(x=axis1, y=axis2, 
                 colour=Mixing_Oxygen), 
             size=5) +
  BioR.theme +
  ggsci::scale_colour_npg() +
  coord_fixed(ratio=1)+
  scale_color_manual(values=c("orange", "blue","green","red"))+
  geom_label_repel(data = sites.long1, 
             aes(label= Date.YYYY.MM.DD, x=axis1, y=axis2),
             alpha=0.7)

plotgg1

# This is Figure S3
ggsave("Figures/NMDS_MAGs.pdf", width = 8, height = 8)

## PHAGE NMDS
virus.taxo <- readRDS("Data/virusTaxo.rds")
virus.taxo1 <- virus.taxo  %>% mutate(Taxonomy = paste(order, family, subfamily)) %>%
  group_by(Taxonomy, Date, Depth, Type) %>% tally()

library(lubridate)
virus.taxo1$Date <- mdy(virus.taxo1$Date)
virus.taxo1$Depth <- str_replace(virus.taxo1$Depth, ".0","")
virus.taxo1$Depth <- paste0(virus.taxo1$Depth, "m")
virus.taxo1$Sample <- paste0(virus.taxo1$Date,"_",
                            virus.taxo1$Depth,"_", virus.taxo1$Type)

virus.taxo1 <- virus.taxo1 %>% select(-Date, -Depth, -Type) %>%
  spread(Sample, n)
virus.taxo1 <- virus.taxo1[,-1:-2]



virus.taxo1 <- virus.taxo1 %>%
  group_by(Taxonomy) %>%
  summarise_if(
    is.numeric,
    sum,
    na.rm = TRUE
  )

rownames.virus.taxo1 <- virus.taxo1$Taxonomy
virus.taxo1 <- virus.taxo1 %>% select(-Taxonomy)
row.names(virus.taxo1) <- rownames.virus.taxo1

virus.taxo1 <- t(virus.taxo1)


Ordination.model2 <- metaMDS(virus.taxo1, distance='bray', 
                             k=2, trymax=1, 
                             autotransform=TRUE,
                             noshare=0.1, expand=TRUE, trace=1, plot=FALSE)

plot2 <- ordiplot(Ordination.model2, choices=c(1,2))

sp.env1 <- sp.env %>% mutate(Sample = paste0(row.names(sp.env),"_Bacteria"),
                             SampleType = "Metagenome")

sp.env2 <- sp.env %>% mutate(Sample = paste0(row.names(sp.env),"_Virus"),
                             SampleType = "Virome")

sp.env.phage <- rbind(sp.env1, sp.env2)


rownames(sp.env.phage) <- sp.env.phage$Sample

row.names(virus.taxo1)

sp.env.phage1 <- sp.env.phage %>% filter(Sample != "2020-08-05_15m_Virus") %>%
  filter(Sample != "2020-08-25_15m_Virus")

nrow(sp.env.phage1)

sites.long2 <- sites.long(plot2, 
                          env.data=sp.env.phage1)

head(sites.long2)

sites.long2$Pairs = paste(sites.long2$Date.YYYY.MM.DD, sites.long2$Sample.Depth.meters)

plotgg2 <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab("NMDS1") +
  ylab("NMDS2") +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
  geom_point(data=sites.long2, 
             aes(x=axis1, y=axis2, color=Mixing_Oxygen, 
                 shape=SampleType), 
             size=5) +
  BioR.theme +
  ggsci::scale_colour_npg() +
  coord_fixed(ratio=1)+
  scale_color_manual(values=c("orange", "blue","green","red"))+
  scale_shape_manual(values=c(16,1))

plotgg2

plotgg3 <- ggplot() + 
  geom_vline(xintercept = c(0), color = "grey70", linetype = 2) +
  geom_hline(yintercept = c(0), color = "grey70", linetype = 2) +  
  xlab("NMDS1") +
  ylab("NMDS2") +  
  scale_x_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +
  scale_y_continuous(sec.axis = dup_axis(labels=NULL, name=NULL)) +    
  geom_point(data=sites.long2, 
             aes(x=axis1, y=axis2, color=Pairs, 
                 shape=SampleType), 
             size=5) +
  BioR.theme +
  #ggsci::scale_colour_npg() +
  coord_fixed(ratio=1)+
  #scale_color_manual(values=c("orange", "blue","green","red"))+
  scale_shape_manual(values=c(16,1))



plotgg3

# This is Figure S3
ggsave(plot  = plotgg2,
       filename = "Figures/NMDS_Phages.pdf", width = 8, height = 8)

library(ggpubr)

ggarrange(plotgg1, plotgg2, labels=c("A","B"),
          nrow=2)

ggsave("Figures/FigureS3_NMDS_MAGs_phages.panel.pdf", width = 8, height = 10)
