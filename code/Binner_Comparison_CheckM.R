# Binner comparison

binner_compare <- read.csv("Data/binner_comparison.csv")
binner_compare <- binner_compare %>% filter(!grepl("renamed", Binner)) %>%
  filter(!is.na(Max))

library(tidyverse)

ggplot(binner_compare, aes(x=Method, y=Average, col=Method))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x =element_text(angle=90, hjust=1),
        strip.background = element_blank(),
        panel.grid = element_blank())+
  facet_wrap(Binner~.)+
  ylab("Average genome quality\nper sample (% completeness)")+
  ggtitle("Comparison of cross-mapping methods with 3 binners ran through metawrap, genome qual by checkM",
          subtitle="Freshwater lake metagenomes, Patricia Tran, Sept 2022")

ggsave("Figures/Comparison_Binners_Genome_qual.pdf", width = 11, height = 6)


ggplot(binner_compare, aes(x=Binner, y=Average, col=Method))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x =element_text(angle=90, hjust=1),
        strip.background = element_blank(),
        panel.grid = element_blank())+
  facet_wrap(Layer~.)+
  xlab("Method and binner")+
  ylab("Average genome quality\nper sample (% completeness)")+
  ggtitle("Comparison of cross-mapping methods with 3 binners ran through metawrap, genome qual by checkM",
          subtitle="Freshwater lake metagenomes")+
  geom_hline(yintercept=50, color="grey", linetype="dotted")



ggsave("Figures/Comparison_Binners_Genome_qual.2.png", 
       width = 11, height = 6)

binner_compare %>% group_by(Layer) %>% tally()
binner_compare %>% group_by(Method) %>% tally() %>% as_tibble()


# Minimum:

ggplot(binner_compare, aes(x=Method, y=Min))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x =element_text(angle=90))+
  facet_wrap(Binner~.)

#Max
ggplot(binner_compare, aes(x=Method, y=Max))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x =element_text(angle=90))+
  facet_wrap(Binner~.)

ggplot(binner_compare, aes(x=Binner, y=Max, col=Method))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x =element_text(angle=90, hjust=1),
        strip.background = element_blank(),
        panel.grid = element_blank())+
  facet_wrap(Layer~.)+
  xlab("Method and binner")+
  ylab("Maximum genome quality\nper sample (% completeness)")+
  ggtitle("Comparison of cross-mapping methods with 3 binners ran through metawrap, genome qual by checkM",
          subtitle="Freshwater lake metagenomes")
  #geom_hline(yintercept=50, color="grey", linetype="dotted")

