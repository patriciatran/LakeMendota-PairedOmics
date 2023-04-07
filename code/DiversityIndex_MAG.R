library(tidyverse)
library(vegan)
library(plyr)

mag_abundance_matrix <- read.csv("Data/MAG_abundance_MATRIX_scaffold_summary_mapping.join_MAG.csv") %>% filter(Final_Bin_Set == "yes")

mag_abundance_matrix.small <- mag_abundance_matrix %>% select(NormCov, Sample.y, Class)

mag_abundance_matrix.small_by_Class <- mag_abundance_matrix.small %>% 
  group_by(Class, Sample.y) %>%
  mutate(abundance_of_class =sum(NormCov))

mag_abundance_matrix.small_by_Class <- mag_abundance_matrix.small_by_Class %>% select(-NormCov) %>% distinct()

sp.matrix <- mag_abundance_matrix.small_by_Class %>% spread(key=Class, value=abundance_of_class)

data(BCI)

H <- diversity(sp.matrix[,-1])
H
simp <- diversity(sp.matrix[-1], "simpson")
simp
shan <- diversity(sp.matrix[-1], "shannon")
shan

pairs(cbind(H, simp, shan), pch="+", col="blue")


### https://www.flutterbys.com.au/stats/tut/tut13.2.html
# Species richness


apply(sp.matrix[,-1]>0,1,sum)

colnames(sp.matrix)[1] <- "Sites"




# Classes  Richness
ddply(sp.matrix,~Sites,function(x) {
  data.frame(RICHNESS=sum(x[-1]>0))
  })

## IN OXIC LAYER:
sample.desc <- read.csv("Data/MT_GaID_metadata.csv") %>% 
  select(Oxygen, Metagenome.sample.name)

sample.desc

# But with microbiome data there are low abundance stuff
richness.plot <- ddply(sp.matrix,~Sites,function(x) {
  data.frame(RICHNESS=sum(x[-1]>0))
}) %>% 
  left_join(x=., y=sample.desc, by=c("Sites"="Metagenome.sample.name"))%>%
  ggplot(aes(x=Sites, y=RICHNESS, col=Oxygen)) +
  geom_point()+
  geom_line(aes(group=Oxygen))+
  geom_vline(xintercept =13.5)+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text= element_text(size=15, color ="black"),
        axis.text = element_text(color="black"),
        axis.title = element_text(face="bold"),
        panel.border = element_rect(colour = "black",linetype = 1,size = 1.5),
        axis.ticks = element_line(size=1.5),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )+
  ylab("Richness")

richness.plot
# ABUNDANCE:

sp.matrix$Sites <- str_replace(sp.matrix$Sites,"ME_","")
sp.matrix$Sites <- str_replace(sp.matrix$Sites,"_Bacteria.*","")
sp.matrix


sample.metadata <- read.csv("Data/Samples-mendota.csv")

cols_oxygen <- c("oxic"="#48AA72",
                 "oxycline"="red",
                 "anoxic"="#3E1E53")

abundance.plot <- ddply(sp.matrix,~Sites,function(x) {
  data.frame(ABUNDANCE=sum(x[-1]))
  }) %>% 
  left_join(x=., y=sample.metadata, by=c("Sites"="Sample.Name"))%>%
  ggplot(aes(x=ymd(Date.YYYY.MM.DD), y=ABUNDANCE, col=Oxygen.level, 
             shape=factor(Sample.Depth.meters))) +
  geom_point()+
  scale_color_manual(values = cols_oxygen)+
  #geom_line(aes(group=Oxygen.level))+
  #facet_grid(Sample.Depth.meters~., scales = "free", space = "free")+
  geom_vline(xintercept =13.5)+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text= element_text(size=12, color ="black"),
        axis.text = element_text(color="black"),
        #axis.title = element_text(face="bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )+
  ylab("Abundance")+
  scale_x_date(breaks = unique(ymd(sample.metadata$Date.YYYY.MM.DD)),
               date_labels = "%b %d")
  
  


abundance.plot 

shannon.plot <- ddply(sp.matrix,~Sites,function(x) {
 data.frame(SHANNON=diversity(x[-1], index="shannon"))
})%>%
  left_join(x=., y=sample.metadata, by=c("Sites"="Sample.Name"))%>%
  ggplot(aes(x=ymd(Date.YYYY.MM.DD), y=SHANNON, col=Oxygen.level, 
             shape=factor(Sample.Depth.meters))) +
  geom_point(size=3)+
  scale_color_manual(values = cols_oxygen)+
  geom_line()+
  #facet_grid(Sample.Depth.meters~., scales = "free", space = "free")+
  #geom_vline(xintercept =13.5)+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text= element_text(size=12, color ="black"),
        axis.text = element_text(color="black"),
        #axis.title = element_text(face="bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1),
        legend.position= "bottom"
  )+
  scale_x_date(breaks = unique(ymd(sample.metadata$Date.YYYY.MM.DD)),
               date_labels = "%b %d")+
  ylab("Shannon-Wiener Index (H')")+
  xlab("Sampling date")

shannon.plot

ggsave("Figures/ShannonIndex.pdf", width = 8, height = 6)

#Shannon-Wiener Index (H′)

shannon.to.plot <- ddply(sp.matrix,~Sites,function(x) {
  data.frame(TRUE_SHANNON=exp(diversity(x[-1], index="shannon")))
})%>%
  left_join(x=., y=sample.metadata, by=c("Sites"="Sample.Name"))

#shannon.to.plot$Sites <- str_replace(shannon.to.plot$Sites, "ME_","")
#shannon.to.plot$Sites <- str_replace(shannon.to.plot$Sites, "_Bacteria_A","")


true_shannon.plot <- shannon.to.plot %>%
  #ggplot(aes(x=Sites, y=TRUE_SHANNON, fill=Oxygen.level)) +
  ggplot(aes(x=ymd(Date.YYYY.MM.DD), y=TRUE_SHANNON, col=Oxygen.level, 
             shape=factor(Sample.Depth.meters))) +
  geom_point()+
  scale_color_manual(values = cols_oxygen)+
  #geom_line(aes(group=Oxygen.level))+
  #facet_grid(Sample.Depth.meters~., scales = "free", space = "free")+
  geom_vline(xintercept =13.5)+
  theme_bw()+
  theme(panel.background = element_blank(),
        panel.grid = element_blank(),
        text= element_text(size=12, color ="black"),
        axis.text = element_text(color="black"),
        #axis.title = element_text(face="bold"),
        axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)
  )+
  scale_x_date(breaks = unique(ymd(sample.metadata$Date.YYYY.MM.DD)),
               date_labels = "%b %d")+
  ylab("True Shannon-Wiener Index")


true_shannon.plot


ggarrange(shannon.plot, true_shannon.plot, common.legend = TRUE,
          nrow=2)

ggsave("Figures/Shannon_diversity_bacterial.Overtime.pdf", width = 10, height = 5)




