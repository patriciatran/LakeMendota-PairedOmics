# Data exploration of checkm results all  MAGs

checkm <- read.csv("Data/checkm_on_refined_bins.csv")

library(tidyverse)

checkm %>% group_by(Marker.lineage) %>% tally()

ggplot(checkm, aes(x=Completeness, y=Contamination))+
  geom_point()+
  geom_vline(xintercept = 50)+
  geom_vline(xintercept = 90, color="red")+
  geom_hline(yintercept=10)+
  geom_hline(yintercept = 5, color="red")+
  theme_bw()+
  theme(text = element_text(size=20))

ggsave("Figures/checkm_1145mags.pdf", width = 10, height = 10)
