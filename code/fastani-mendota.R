#  Comparing ANI between the MAGs in this dataset and previously published MAGs from Lake Mendota:
# Patricia : Patricia Tran
# BP: Ben  Peterson
# EAM: Elizabeth McDaniel
# AL: Alex Linz


library(tidyverse)

mendota.compare <- read.csv("Data/fastANI_other_Mendota-mags.csv", header=TRUE)

mendota.compare.table <- mendota.compare %>% group_by(MAG.ID.Patricia, Phylum, Class, Order, Previous.Ref) %>% 
  tally() %>%
  spread(Previous.Ref, n)

mendota.compare.table <- mendota.compare.table %>% mutate(sumMAG = sum(AL,BP,EAM, na.rm=TRUE))

mendota.compare %>% group_by(Phylum) %>% tally()

ggplot(mendota.compare, aes(x=Order, y=ANI.above.80., col=Previous.Ref))+
  geom_boxplot()+
  theme_bw()+
  theme(axis.text.x = element_text(angle=90, color="black"))+
  scale_color_manual(values=c("red","blue","dark green"))+
  xlab("Taxonomy of MAG in Tran et al. 202X")+
  ylab("% ANI similarity to previous dataset")+
  geom_hline(yintercept = 97)+
  annotate("text",label="97% line", x=1, y=101)

ggsave("Figures/FastANI_Patricia_vs_previousstudentsMAGS_Class.pdf", width = 11, height = 8.5)

write.table(mendota.compare.table, "Data/mendota.compare.table.fastani.tsv", sep="\t", row.names = FALSE, quote = FALSE)

mendota.compare  %>% distinct(MAG.ID.Patricia) %>% nrow()
mendota.compare %>% filter(ANI.above.80. > 80) %>% distinct(MAG.ID.Patricia) %>% nrow()
mendota.compare %>% filter(ANI.above.80. > 97) %>% distinct(MAG.ID.Patricia) %>% nrow()
mendota.compare %>% filter(ANI.above.80. > 98) %>% distinct(MAG.ID.Patricia) %>% nrow()
mendota.compare %>% filter(ANI.above.80. > 99) %>% distinct(MAG.ID.Patricia) %>% nrow()
mendota.compare %>% filter(ANI.above.80. > 99.5) %>% distinct(MAG.ID.Patricia) %>% nrow()
test <- mendota.compare %>% filter(ANI.above.80. > 99.8)

test
