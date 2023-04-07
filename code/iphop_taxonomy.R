iphop.results <- read.csv("Data/Host_prediction_to_genome_m90.csv")


library(tidyverse)

iphop.results$Host.taxonomy[1]
str(iphop.results)
iphop.results <- iphop.results %>% separate(col = Host.taxonomy, sep = ";", remove = FALSE,
                          c("domain","phylum","class","order","family","genus","species"))

# iphop.results %>% group_by(phylum, class) %>% tally() %>% arrange(desc(n)) %>% 
#   ggplot(aes(x=n, y=class))+
#   facet_grid(phylum~.)+
#   geom_col()

MAG.final <- read.table("Data/FinalBinSet/MAG_Characteristics.tsv", sep="\t", header=TRUE)

iphop.results <- iphop.results %>% mutate(host.in.MAG = ifelse(Host.taxonomy %in% MAG.final$classification,
                         yes="In binned MAG",
                         no="not in binned MAG"))

iphop.results %>% group_by(host.in.MAG) %>% tally()
iphop.results %>% distinct(Virus)
iphop.results %>% group_by(domain, phylum, class, host.in.MAG, Main.method) %>% 
  tally() %>%
  ggplot(aes(x=n, 
             y=class, 
             fill=host.in.MAG))+
  geom_col()+
  theme_bw()+
  xlab("Number of phages")+
  ylab("host taxonomy class")+
  scale_fill_manual(values=c("black","grey"))+
  facet_grid(domain~Main.method, scales="free",space = "free")+
  theme(panel.grid = element_blank(),
        strip.background = element_blank(),
        strip.text.y = element_text(angle=0),
        strip.text.x = element_text(angle=90),
        text = element_text(color="black", size=20),
        legend.position = "bottom")

ggsave("Figures/iPhop_host_taxonomy_distribution.pdf", width = 11, height = 11)

iphop.results %>% group_by(phylum, class, host.in.MAG, Main.method, Confidence.score) %>% tally() %>%
  ggplot(aes(x=class, y=Confidence.score, col=host.in.MAG))+
  geom_violin()+
  geom_point()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90))+
  scale_color_manual(values=c("red","blue"))






