# Patricia Tran
# Phage host linkages, from the  iPhop results

links <- read.csv("~/Downloads/Host_prediction_to_genome_m90.csv")

library(tidyverse)
unique(links$Host.genome)

links.mags.only <- links %>% filter(grepl("Ga", Host.genome))

how.many.virus.for.each.mag <- links.mags.only %>% group_by(Host.genome) %>% tally()
how.many.mags.for.each.virus <- links.mags.only %>% group_by(Virus) %>% tally()

hist(how.many.mags.for.each.virus$n, breaks= 10)

hist(how.many.virus.for.each.mag$n, breaks = 10)

test <- links.mags.only %>% group_by(Virus, Host.taxonomy) %>% tally()
test

write.table(test, "Data/Phage-specificity.tsv", sep="\t", quote = FALSE, row.names = FALSE)

test2 <- links.mags.only %>% filter(grepl("Ga",Host.genome)) %>% group_by(Virus, Host.taxonomy) %>% tally()


links.mags.only.pasted.viruses <- links.mags.only %>% select(Host.genome, Virus) %>% 
  group_by(Host.genome) %>% mutate(Phages = paste0(Virus, sep=";", collapse = "")) %>% distinct() %>% select(-Virus)

links.mags.only.pasted.viruses <- left_join(links.mags.only.pasted.viruses, how.many.virus.for.each.mag)

links.mags.only.pasted.viruses.with.taxo <- links.mags.only %>% select(Host.genome, Host.taxonomy)%>%
  left_join(links.mags.only.pasted.viruses) %>% distinct()

write.table(links.mags.only.pasted.viruses.with.taxo, "Data/Host-Phage.pairs.taxonomy.tsv", sep="\t", row.names = FALSE, quote = FALSE)


PHAGE.abund <- read.table("Data/to_share/phages_abundance_matrix_metagenomes.tsv", sep="\t",header=TRUE)
MAG.abund <- read.table("Data/to_share/MAG_abundance_matrix_rel_abund.tsv", sep="\t", header=TRUE)
MAG.expression <- read.table("Data/to_share/MAG_RPKM_normalized.tsv", sep="\t", header=TRUE)
Phage.expression <- read.table("Data/to_share/phages_rpkm_normalized_matrix.tsv", sep="\t", header=TRUE)

MAG.abund
MAG.abund <- MAG.abund %>% filter(Genome != "unmapped")

library(ggpubr)
library(lubridate)

sample.metadata <- read.csv("Data/to_share/Samples-mendota.csv")

cols_oxygen <- c("oxic"="#48AA72",
                 "oxycline"="#48658C",
                 "anoxic"="#3E1E53")

shape_manual <- c("MAG"=16,
                  "PHAGE"=4)

View(links.mags.only.pasted.viruses.with.taxo)
#list.of.i <- c(30, 31, 5, 23, 1, 42) # these are the host that have a phage that is expressed multiple times.
list.of.i <- sort(list.of.i)

sample.metadata$Mixing.regime <- factor(sample.metadata$Mixing.regime,
                                        levels=c("stratified","mixed"))

sample.metadata$Oxygen.level <- factor(sample.metadata$Oxygen.level,
                                       levels=c("oxic","oxycline","anoxic"))
samples.ordered <- sample.metadata %>% arrange(Mixing.regime, 
                                               Oxygen.level,
                                               Date.YYYY.MM.DD, 
                                               Sample.Depth.meters) %>% select(Sample.Name)
samples.ordered$Sample.Name

sample.metadata$Sample.Name <- factor(sample.metadata$Sample.Name,
                                      levels=samples.ordered$Sample.Name)

sample.metadata$Sample.Name

levels(sample.metadata$Sample.Name)

for (i in 1:nrow(links.mags.only.pasted.viruses.with.taxo)){
  print(i)
  #taxonomy.to.print <- "Ga0485159_metabat2_jgi.002"
  #mag.to.take <- "Ga0485159_metabat2_jgi.002"
  taxonomy.to.print <- links.mags.only.pasted.viruses.with.taxo$Host.taxonomy[i]
  mag.to.take <- links.mags.only.pasted.viruses.with.taxo$Host.genome[i]
  
  
  
  phage.to.get <- links.mags.only.pasted.viruses.with.taxo$Phages[i] %>% strsplit(";")
  
  phage.to.get <- phage.to.get[[1]]
  phage.to.get
  
  
  mag.abund <- MAG.abund %>% filter(Genome == mag.to.take) %>%
    gather("sample","relabund",2:17) %>% mutate(type = "MAG", 
                                                genome = Genome)
  
  mag.expression <- MAG.expression %>% select(-2,-3) %>% filter(Genome == mag.to.take) %>%
    gather("sample","rpkm_norm",2:17) %>% mutate(type = "MAG", 
                                                 genome = Genome)
  
  
  phage.abund <- PHAGE.abund %>% filter(Phage %in% phage.to.get) %>%
    gather("sample","relabund", 2:17) %>% mutate(type = "PHAGE",
                                                 genome = Phage)
  
  phage.expression <- Phage.expression %>% filter(Phage %in% phage.to.get) %>%
    gather("sample","rpkm_norm", 2:17) %>% mutate(type = "PHAGE",
                                                  genome = Phage)
  ## abundance
  full_join_mag_phage <- full_join(mag.abund, phage.abund) %>% select(-Genome, -Phage)
  
  full_join_mag_phage$sample <- str_replace(full_join_mag_phage$sample, "X","")
  full_join_mag_phage$sample <- str_replace(full_join_mag_phage$sample, "\\.","-")
  full_join_mag_phage$sample <- str_replace(full_join_mag_phage$sample, "\\.","-")
  
  full_join_mag_phage <- left_join(full_join_mag_phage, sample.metadata, by=c("sample"="Sample.Name"))
  full_join_mag_phage$Date.YYYY.MM.DD <- ymd(full_join_mag_phage$Date.YYYY.MM.DD)
  
  full_join_mag_phage$Oxygen.level <- factor(full_join_mag_phage$Oxygen.level,
                                             levels=c("oxic","oxycline","anoxic"))
  
  
  ## expression:
  full_join_mag_phage_expression <- full_join(mag.expression, phage.expression) %>% select(-Genome, -Phage)
  
  full_join_mag_phage_expression$sample <- str_replace(full_join_mag_phage_expression$sample, "X","")
  full_join_mag_phage_expression$sample <- str_replace(full_join_mag_phage_expression$sample, "\\.","-")
  full_join_mag_phage_expression$sample <- str_replace(full_join_mag_phage_expression$sample, "\\.","-")
  
  full_join_mag_phage_expression <- left_join(full_join_mag_phage_expression, sample.metadata, by=c("sample"="Sample.Name"))
  full_join_mag_phage_expression$Date.YYYY.MM.DD <- ymd(full_join_mag_phage_expression$Date.YYYY.MM.DD)
  
  full_join_mag_phage_expression$Oxygen.level <- factor(full_join_mag_phage_expression$Oxygen.level,
                                                        levels=c("oxic","oxycline","anoxic"))
  
  
  
  
  panelA <- ggplot(full_join_mag_phage %>% filter(type=="MAG"), aes(x=Date.YYYY.MM.DD, 
                                                                    y=relabund, 
                                                                    group=type, 
                                                                    col=Oxygen.level
                                                                    ))+
    geom_point(size=2)+
    geom_line(alpha=1, size=0.5, col="black")+
    facet_grid(Sample.Depth.meters~., scales = "free_y")+
    theme_bw()+
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1),
          panel.grid.minor = element_blank(),
          plot.subtitle=element_text(size=8, hjust=0, 
                                     face="italic", 
                                     color="black"),
          strip.text.y = element_text(angle=0))+
    ylab("Relative abundance amongs all reads\n (inculding unmmapped)")+
    xlab("Sampling date")+
    ggtitle("Host abundance")+
    scale_color_manual(values= cols_oxygen)+
    scale_shape_manual(values=shape_manual)+
    scale_x_date(breaks = unique(ymd(sample.metadata$Date.YYYY.MM.DD)))
  
  panelA
  

  
  panelB <- ggplot(full_join_mag_phage %>% filter(type == "PHAGE") %>%
                     group_by(Date.YYYY.MM.DD, Sample.Depth.meters) %>%
                     mutate(Mean = mean(relabund),
                            Sd = sd(relabund)),
                   aes(x=Date.YYYY.MM.DD, 
                       y=Mean, 
                       group=type, 
                       col=Oxygen.level,
                   ))+
    geom_errorbar(aes(ymin=Mean-Sd,
                      ymax=Mean+Sd))+
    geom_point(size=2)+
    #geom_line(alpha=0.5, size=0.5, col="black")+
    facet_grid(Sample.Depth.meters~., scales = "free_y")+
    theme_bw()+
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1),
          panel.grid.minor = element_blank(),
          plot.subtitle=element_text(size=8, hjust=0, 
                                     face="italic", 
                                     color="black"),
          strip.text.y = element_text(angle=0))+
    ylab("Relative abundance amongs all reads\n (inculding unmmapped)")+
    xlab("Sampling date")+
    ggtitle("Phage abundance")+
    scale_color_manual(values= cols_oxygen)+
    #scale_shape_manual(values=shape_manual)+
    scale_x_date(breaks = unique(ymd(sample.metadata$Date.YYYY.MM.DD)))
  
  panelB 
  
  panelC <-ggplot(full_join_mag_phage_expression %>% filter(type == "MAG"), aes(x=Date.YYYY.MM.DD, 
                                                                                y=rpkm_norm, 
                                                                                group=type, 
                                                                                col=Oxygen.level,
                                                                                shape=type))+
    geom_point(size=2)+
    geom_line(alpha=1, size=0.5, col="black")+
    facet_grid(Sample.Depth.meters~., scales = "free_y")+
    theme_bw()+
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1),
          panel.grid.minor = element_blank(),
          plot.subtitle=element_text(size=8, hjust=0, face="italic", color="black"),
          strip.text.y = element_text(angle=0))+
    ylab("RPKM normalized by internal standard")+
    xlab("Sampling date")+
    ggtitle("Host expression")+
    scale_color_manual(values= cols_oxygen)+
    scale_shape_manual(values=shape_manual)+
    scale_x_date(breaks = unique(ymd(sample.metadata$Date.YYYY.MM.DD)))
  
  panelC
  
  panelD <- ggplot(full_join_mag_phage_expression %>% filter(type == "PHAGE") %>%
                     group_by(Date.YYYY.MM.DD, Sample.Depth.meters) %>%
                     mutate(Mean = mean(rpkm_norm),
                            Sd = sd(rpkm_norm)),
                   aes(x=Date.YYYY.MM.DD, 
                       y=Mean, 
                       group=type, 
                       col=Oxygen.level,
                   ))+
    geom_errorbar(aes(ymin=Mean-Sd,
                      ymax=Mean+Sd))+
    geom_point(size=2)+
    #geom_line(alpha=0.5, size=0.5, col="black")+
    facet_grid(Sample.Depth.meters~., scales = "free_y")+
    theme_bw()+
    theme(strip.background = element_blank(),
          axis.text.x = element_text(angle=90, hjust=1),
          panel.grid.minor = element_blank(),
          plot.subtitle=element_text(size=8, hjust=0, face="italic", color="black"),
          strip.text.y = element_text(angle=0))+
    ylab("RPKM normalized by internal standard")+
    xlab("Sampling date")+
    ggtitle("Phage expression")+
    scale_color_manual(values= cols_oxygen)+
    scale_shape_manual(values=shape_manual)+
    scale_x_date(breaks = unique(ymd(sample.metadata$Date.YYYY.MM.DD)))
  
  panelD
  
  if(length(phage.to.get) == 1){
    panelB <- panelB + 
      geom_line(alpha=1, size=0.5, col="black")
    panelB
    
    panelD <- panelD + 
      geom_line(alpha=1, size=0.5, col="black")
    panelD
  }
  else{
    print("no modif")
  }

  
  combined.figure <-ggarrange(panelA, panelB, panelC, panelD, 
                              labels=c("A","B","C","D"),
                              nrow=2, ncol=2, common.legend = TRUE)
  
  
  annotate_figure(combined.figure, 
                  top = text_grob(paste0(mag.to.take, "\n",taxonomy.to.print), 
                                  color = "black", 
                                  size = 12))
  
  
  ggsave(paste0("Figures/Option3_phage-host.patterns_",mag.to.take,".pdf"), 
         width = 7, height = 8.5)
  

  if (i %in% list.of.i){
    # this is one of the MAGs that we want the phage for
    print("one of the phages of interest")
    
    
    full_join_mag_phage.save <- full_join_mag_phage %>% select(sample, relabund, genome) %>%
      spread(genome, relabund) %>%
      arrange(factor(sample, levels = samples.ordered$Sample.Name))
    
    write.table(full_join_mag_phage.save, 
              file=paste0("Data/Full_join_mag_phage_",i,".tsv"),
              sep="\t",
              quote = FALSE,
              row.names = FALSE
              )
    
    full_join_mag_phage_expression.save <- full_join_mag_phage_expression %>% select(sample, rpkm_norm, genome) %>%
      spread(genome, rpkm_norm) %>%
      arrange(factor(sample, levels = samples.ordered$Sample.Name))
    
    full_join_mag_phage_expression.save
    
    write.table(full_join_mag_phage_expression.save, 
                file=paste0("Data/Full_join_mag_phage_expression",i,".tsv"),
                sep="\t",
                quote = FALSE,
                row.names = FALSE
    )
    
  }
  else{print("continue")}
  
}

library(plotrix)
twoord.plot()


         
head(links.mags.only.pasted.viruses.with.taxo)

str_replace(links.mags.only.pasted.viruses.with.taxo$Phages, "_.*;","")
strsplit(links.mags.only.pasted.viruses.with.taxo$Phages, ";")



links.mags.only1 <-links.mags.only %>% select(Host.genome, Virus) 
links.mags.only1$Virus_sample <- str_replace(links.mags.only1$Virus, "_.*","")

samples.metadata <- read.csv("Data/GaID_match.csv")

links.mags.only1$Host.genome_sample <- str_replace(links.mags.only1$Host.genome, "_.*","")

links.mags.only1 <- left_join(links.mags.only1, samples.metadata, by=c("Virus_sample"="GaID"))

library(lubridate)

links.mags.only1$Date <- mdy(links.mags.only1$Date)

unique.genomes <- unique(links.mags.only1$Host.genome)
xaxis <- mdy(unique(samples.metadata$Date))

ggplot(links.mags.only1, aes(x=Date, y=TYPE))+
  geom_tile()+
  facet_wrap(Host.genome~.)

links.mags.only2 <- left_join(links.mags.only1, samples.metadata, by=c("Host.genome_sample"="GaID"))

links.mags.only2$Date.y <- mdy(links.mags.only2$Date.y)

toplot1 <- links.mags.only2 %>% select(Host.genome,Date.x, TYPE.x, Depth.x)

toplot2 <- links.mags.only2 %>% select(Host.genome, Date.y, Depth.y)
toplot2$TYPE.x <- "Representative genome"

colnames(toplot1) 

colnames(toplot2) <- c("Host.genome","Date.x","Depth.x","TYPE.x")
toplot.combined <- full_join(toplot1, toplot2)

toplot.combined$TYPE.x <- factor(toplot.combined$TYPE.x, 
                                 levels=c("Metagenomes","Viromes","Representative genome"))

head(toplot.combined)

toplot.combined <- toplot.combined %>% group_by(Host.genome, Date.x, TYPE.x, Depth.x) %>% tally() %>% distinct()
toplot.combined <- toplot.combined %>% mutate(Count= ifelse(TYPE.x == "Representative genome", 1, n))


head(toplot.combined)

ggplot(toplot.combined, aes(x=Date.x, y=TYPE.x, fill=TYPE.x, 
                            size=Count, col=TYPE.x))+
  geom_point()+
  geom_vline(xintercept = mdy("10-18-2020"), linetype="dashed")+
  facet_wrap(Host.genome~.)+
  theme_bw()+
  theme(panel.grid.minor = element_blank(),
        strip.background = element_blank())+
  scale_fill_manual(values=c("blue","black","black"))

ggsave("Figures/Infection_patterns.pdf", width = 11, height = 7)

host.genomes.unique <- unique(toplot.combined$Host.genome)

levels(toplot.combined$TYPE.x)
for (i in 1:length(host.genomes.unique)){
  print(i)
  
  toplot.combined %>% filter(Host.genome == host.genomes.unique[i]) %>%
    ggplot(aes(x=Date.x, y=-Depth.x, shape=TYPE.x, size=Count))+
    geom_point()+
    scale_shape_manual(values=c(0,2,4),breaks = levels(toplot.combined$TYPE.x))+
    ylim(-23.5,0)+
    theme_bw()+
    theme(axis.text.x = element_text(angle=90),
          panel.grid.minor = element_blank())+
    ylab("Depth (m)")+
    xlab("Sampling date")+
    #scale_x_date(date_breaks = mdy(samples.metadata$Date))+
    xlim(c(min(mdy(samples.metadata$Date)), max(mdy(samples.metadata$Date))))+
    geom_vline(xintercept = mdy("10-18-2020"), linetype = "dashed")+
    ggtitle(host.genomes.unique[i])
  
  ggsave(paste0("Figures/",host.genomes.unique[i],"_phage_host.pdf"), width = 10, height = 5)
}


## Plotting more pannels:
full_join_mag_phage_expression %>% filter(Date.YYYY.MM.DD == "2020-07-24" & Sample.Depth.meters=="23.5")
