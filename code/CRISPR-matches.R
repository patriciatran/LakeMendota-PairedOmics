# some figures for my ISME poster

library(gggenes)
library(ggpubr)
library(tidyverse)
library(lubridate)


crispr.match <- read.csv("Data/CRISPR-spacers-matches-all-phages.csv", header=TRUE)

crispr.match %>% select(MAG.name, MAG.phylum, How.many.different.phages.is.this.MAG.protected.against.) %>%
  distinct() %>%
  ggplot(aes(x=MAG.name,
             y=How.many.different.phages.is.this.MAG.protected.against.,
             fill=MAG.phylum))+
  geom_col()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, hjust=1))+
  ylab("Number of distinct phage \n with perfect CRISPR spacer match")+
  xlab("MAG name")

test <- crispr.match %>% select(Phage.match.in.all.the.VIBRANT.results, How.many.different.phyla.does.this.phage.infect, How.many.different.phages.is.this.MAG.protected.against.) %>%
  distinct()  %>% gather("Question","Value",2:3) 


test %>% 
  ggplot(aes(x=Phage.match.in.all.the.VIBRANT.results,
             y=Value,
             fill=Question))+
  geom_col(position="dodge")


pilercr <- read.csv("Data/Ga0485162_0000441_pilercr_sequences.csv")


str(pilercr)

# Plot 1 for the MAG:

plot2 <- ggplot(pilercr %>% filter(Genome == "Ga0485162_0000441 pilercr annotation")) + 
  geom_gene_arrow(aes(xmin = Start.Coord, 
                      xmax = End.Coord, 
                      y=Genome, 
                      fill=Gene_Name.2)) +
  facet_wrap(.~ Genome, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3")+
  theme_genes()+
  ylab("Phage scaffold name")+
  scale_fill_manual(values=c("red","white","blue","green"),
                    breaks=c("CRISPR Cas","Other","Repeat","Spacer"))+
  ggtitle("Bacterial host, zoom in")

plot2 
#ggsave("Figures/Ga0485162_0000441_pilercr_sequences.pdf", width = 9, height = 2.5)

pilercr$Genome2 <- "Bacteria Ga0485162_0000441"


plot1 <- ggplot(pilercr) +
  geom_gene_arrow(aes(xmin = Start.Coord, 
                      xmax = End.Coord, 
                      y=Genome2, 
                      fill=Gene_Name.2)) +
  #facet_wrap(.~ Genome, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3")+
  theme_genes()+
  ylab("Phage scaffold name")+
  scale_fill_manual(values=c("red","white","blue","green"),
                    breaks=c("CRISPR Cas","Other","Repeat","Spacer"))+
  ggtitle("Bacterial host, full scaffold")
plot1

# Plot 2 for the phage:
phage <- read.table("Data/Phage_Ga0485157_0006781.txt", header=TRUE, sep="\t")

phage$Genome2 <- "Phage"

plot3 <- ggplot(phage) + 
  geom_gene_arrow(aes(xmin = Start.Coord, 
                      xmax = End.Coord, 
                      y=Genome2, 
                      fill=Gene.Name)) +
  #facet_wrap(.~ Genome, scales = "free", ncol = 1) +
  #scale_fill_brewer(palette = "Set3")+
  theme_genes()+
  ylab("Phage scaffold name")+
  scale_fill_manual(values=c("black","white","purple"),
                    breaks=c("100% match to host","Other","Phage protein"))+
  ggtitle("Phage matching the host")

plot3



ggarrange(plot1,plot2,plot3, ncol=1,
          labels=c("A","B","C"))

ggsave("Figures/Phage-host-crispr-spacers-Ga0485162_0000441.pdf", width = 9, height = 5)

# Abundance of the 10 phages:
abun.matrix <- read.csv("Data/MAG_abundance_MATRIX_scaffold_summary_mapping.join_MAG.csv")

list.of.phages.with.matching.crispr.spacers
list.of.mags.with.crispr <- c("Ga0485162_metabat2_ours.026","Ga0485158_metabat2_ours.011_sub","Ga0485169_metabat2_ours.021_sub","Ga0485171_metabat2_ours.182","Ga0485159_metabat2_jgi.002","Ga0485161_metabat2_jgi.003_sub","Ga0485158_metabat2_ours.184","Ga0485169_maxbin.156","Ga0485172_maxbin.002","Ga0485163_metabat1.110_sub")


abund.matrix.subset <- abun.matrix %>% filter(bin %in% list.of.mags.with.crispr)
abund.matrix.subset <- abund.matrix.subset %>% separate(Sample.y, sep="_", into=c("Lake","Date","Depth","Type","Replicate"), remove=FALSE)


ggplot(abund.matrix.subset, aes(x=ymd(Date), y= NormCov, fill=Class))+
  #geom_point(size=2, alpha=0.8)+
  #geom_line()+
  theme_bw()+
  geom_col()+
  facet_grid(as.numeric(str_replace(Depth,"m",""))~.)+
  #theme(legend.position = "bottom")+
  #geom_vline(xintercept=ymd("2020-10-18"))+
  scale_fill_manual(values = c("#3F85C6", #choroflexi
                               "#56B2A8", #bacteroidetes
                               "#8BAE7B", #cyano
                               "#3B2E72", # gamma
                               "#9E4293", # kiri
                               "#A07D46", #bdello
                               "#B486B8"# verruco
                               ))+
  xlab("Date")+
  ylab("Normalized coverage per 100M reads")+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle=0))

ggsave("Figures/Abundance_MAGS_with_CRISPR_matches.pdf", width = 11, height = 8.5)

unique(abund.matrix.subset$Class)

abund.matrix.subset.with.phage <- abun.matrix %>% filter(Class %in% abund.matrix.subset$Class) %>%
  filter(Final_Bin_Set == "yes")

abund.matrix.subset.with.phage <- abund.matrix.subset.with.phage %>%
  mutate(Phage_status = ifelse(bin %in% list.of.mags.with.crispr,
                               "Phage pair found", #if true
                               "Phage pair not found or unclear")) # if false

abund.matrix.subset.with.phage <- abund.matrix.subset.with.phage %>% separate(Sample.y, sep="_", into=c("Lake","Date","Depth","Type","Replicate"), remove=FALSE)

for (i in 1:length(unique(abund.matrix.subset.with.phage$Class))){
  ggplot(abund.matrix.subset.with.phage %>% filter(Class == unique(abund.matrix.subset.with.phage$Class)[i]), 
         aes(x=ymd(Date), y= NormCov, fill=Phage_status))+
    #geom_point(size=2, alpha=0.8)+
    #geom_line()+
    theme_bw()+
    geom_col(position="dodge")+
    facet_grid(as.numeric(str_replace(Depth,"m",""))~.)+
    #theme(legend.position = "bottom")+
    #geom_vline(xintercept=ymd("2020-10-18"))+
    # scale_fill_manual(values = c("#3F85C6", #choroflexi
    #                              "#56B2A8", #bacteroidetes
    #                              "#8BAE7B", #cyano
    #                              "#3B2E72", # gamma
    #                              "#9E4293", # kiri
    #                              "#A07D46", #bdello
    #                              "#B486B8"# verruco
    # ))+
    scale_fill_manual(values=c("black","grey"),
                      breaks=c("Phage pair found","Phage pair not found or unclear"))+
    xlab("Date")+
    ylab("Normalized coverage per 100M reads")+
    theme(panel.grid = element_blank(),
          strip.background = element_blank(),
          strip.text.y = element_text(angle=0),
          legend.position="bottom",
          text=element_text(size=30))+
    ggtitle(unique(abund.matrix.subset.with.phage$Class)[i])
  
  ggsave(paste0("Figures/Classes-mags-with-phages_",unique(abund.matrix.subset.with.phage$Class)[i],".pdf"),
         width = 12, height = 8)
  print("next")
  

}

ggplot(abund.matrix.subset.with.phage, 
       aes(x=ymd(Date), y= NormCov, fill=Phage_status))+
  #geom_point(size=2, alpha=0.8)+
  #geom_line()+
  theme_bw()+
  geom_col(position="dodge")+
  facet_grid(as.numeric(str_replace(Depth,"m",""))~Class)+
  #theme(legend.position = "bottom")+
  #geom_vline(xintercept=ymd("2020-10-18"))+
  # scale_fill_manual(values = c("#3F85C6", #choroflexi
  #                              "#56B2A8", #bacteroidetes
  #                              "#8BAE7B", #cyano
  #                              "#3B2E72", # gamma
  #                              "#9E4293", # kiri
  #                              "#A07D46", #bdello
  #                              "#B486B8"# verruco
  # ))+
  scale_fill_manual(values=c("black","grey"),
                    breaks=c("Phage pair found","Phage pair not found or unclear"))+
  xlab("Date")+
  ylab("Normalized coverage per 100M reads")+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle=0),
        legend.position="bottom",
        text=element_text(size=12))

ggsave("Figures/mags_abundance_phage_vs_no_phages.pdf", width = 12, height = 7)



