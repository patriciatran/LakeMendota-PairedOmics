compare_binning <- read.csv("Data/Comparingbinningmethods.csv")

library(tidyverse)

ggplot(compare_binning, aes(x=Sample, y=Number.of.bins, fill=Method))+
  geom_col(position="dodge")+
  #facet_wrap(Method~.)+
  theme_bw()+
  xlab("Number of MAGs obtained")+
  theme(axis.text.x = element_text(angle=90),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        text=element_text(size=12, color="black"))+
  scale_fill_manual(values=c("#55D6BE","#7D5BA6","#FC6471"))
  
ggsave("Figures/Comparing_binning_methods.pdf", width = 12, height = 7)

ggplot(compare_binning)+
  geom_point(aes(x=Software, y=Number.of.bins, col=Sample))+
  geom_path(aes(x=Software, y=Number.of.bins, col=Sample, group=Sample))+
  facet_wrap(.~Method)+
  theme_bw()+
  xlab("Binning method (all through metawrap)")+
  ylab("Number of bins obtained")+
  theme(axis.text.x = element_text(angle=90),
        strip.background = element_blank(),
        panel.grid = element_blank(),
        text=element_text(size=15, color="black"))

ggsave("Figures/Comparing_binning_methods_line_plot.pdf", width = 12, height = 7)


  

