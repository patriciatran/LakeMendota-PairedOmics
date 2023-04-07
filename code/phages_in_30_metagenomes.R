vibrant_results <- read.csv("Data/VIBRANT-results.csv")

library(tidyverse)

my_fave_theme <-   theme(panel.background = element_blank(),
                         panel.grid = element_blank(),
                         text= element_text(size=15, color ="black"),
                         axis.text = element_text(color="black"),
                         axis.title = element_text(face="bold"),
                         #panel.border = element_rect(colour = "black",linetype = 1,size = 1.5),
                         axis.ticks = element_line(size=1.5)
)

dim(vibrant_results)

vibrant_results %>% group_by(Type) %>% tally()


all.phages.plot <- ggplot(vibrant_results, aes(x=Type, 
                                               y=Putative.phages.found))+
  #geom_boxplot()+
  #geom_violin()+
  
  geom_point(alpha=0.5)+
  theme_bw()+
  my_fave_theme +
  ylab("Putative phages found")+
  xlab("Sample Type")+
  stat_summary(fun="mean", color="red")+
  stat_summary(fun="median", color="blue")



all.phages.plot


lytic.plot <- ggplot(vibrant_results, aes(x=Type, y=Lytic))+
  #geom_boxplot()+
  geom_violin()+
  theme_bw()+
  my_fave_theme +
  ylab("Lytic phages")+
  xlab("Sample Type")+
  geom_point(alpha=0.5)+
  stat_summary(fun="mean", color="red")+
  stat_summary(fun="median", color="blue")



lysogenic.plot <- ggplot(vibrant_results, aes(x=Type, y=Lysogenic))+
  #geom_boxplot()+
  geom_violin()+
  theme_bw()+
  my_fave_theme +
  ylab("Lysogenic phages")+
  xlab("Sample Type")+
  geom_point(alpha=0.5)+
  stat_summary(fun="mean", color="red")+
  stat_summary(fun="median", color="blue")


percent.of.data.phages.plot <- ggplot(vibrant_results, aes(x=Type, y=Putative.phages.found.percent))+
  geom_point()+
  geom_violin()+
  theme_bw()+
  my_fave_theme +
  ylab("Percent of sequences in \n metagenome that are putative phages")+
  xlab("Sample Type")+
  scale_y_continuous(labels = scales::percent)+
  stat_summary(fun="mean", color="red")+
  stat_summary(fun="median", color="blue")



# The two groups are dependent on each other: Paired samples t-test
#t.test(x, y, paired = TRUE, alternative = "two.sided")
#install.packages("PairedData")
library(PairedData)


metagenomes<- subset(vibrant_results %>% filter(Paired.Sample.Number != "Unpaired"), 
                     Type == "Metagenome", Putative.phages.found.percent,
                              drop = TRUE)

viromes<- subset(vibrant_results %>% filter(Paired.Sample.Number != "Unpaired"),  Type == "Virome", Putative.phages.found.percent,
                     drop = TRUE)

paired.data <- paired(metagenomes,viromes)
head(paired.data)

library(ggpubr)
plotA <- ggpaired(paired.data, cond1="metagenomes", 
         cond2="viromes")+
  ylab("Percent of data that are phages")+
  xlab("Sample type")

plotA 
  
ggsave(plotA, "Figures/Pairs.phages.boxplot.pdf", units = "in", height = 8, width = 8)

# compute the difference
d <- with(vibrant_results %>% filter(Paired.Sample.Number != "Unpaired"), 
          Putative.phages.found.percent[Type == "Metagenome"] - Putative.phages.found.percent[Type == "Virome"])

d
# Shapiro-Wilk normality test for the differences
shapiro.test(d) # => p-value = p-value = 0.8523



#From the output, the p-value is greater than the significance level 0.05 implying that the distribution of the differences (d) are not significantly different from normal distribution. In other words, we can assume the normality.
res <- t.test(metagenomes, viromes, paired = TRUE)
res

# Paired t-test
# 
# data:  metagenomes and viromes
# t = -11.236, df = 13, p-value = 4.581e-08
# alternative hypothesis: true difference in means is not equal to 0
# 95 percent confidence interval:
#   -0.5896818 -0.3994911
# sample estimates:
#   mean of the differences 
# -0.4945864 

percent.of.data.phages.plot  = ggpaired(paired.data, cond1="metagenomes", 
         cond2="viromes")+
  ylab("% of dataset \nthat are phages")+
  xlab("Sample type")+
  my_fave_theme+
  annotate(geom="text", label="* p=4.581e-08", y=1, x=1.5)+
  annotate(geom="segment", x=1, xend=2, y=0.90, yend=0.90)+
  scale_y_continuous(labels = scales::percent)
  

percent.of.data.phages.plot

# Figure 4A
ggsave("Figures/Pairs.phages.boxplot.pvalue.pdf", units = "in", height = 8, width = 8)


prophages.plot <- ggplot(vibrant_results, 
                aes(x=Type, y=Prophages.count))+
  geom_violin()+
  geom_point(alpha=0.5)+
  theme_bw()+
  my_fave_theme +
  ylab("Prophages")+
  xlab("Sample Pair")+
  stat_summary(fun="mean", color="red")+
  stat_summary(fun="median", color="blue")

library(cowplot)

top_row <- plot_grid(percent.of.data.phages.plot, all.phages.plot, 
                     labels=c('A', 'B'), 
                     label_size=12, ncol=2)
top_row

bottom_row <- plot_grid(lytic.plot, lysogenic.plot, prophages.plot, 
                        labels = c('C','D', 'E'), 
                        label_size = 12, ncol=3)
bottom_row

# checkv_qualityall.R 

plot_grid(top_row, 
          bottom_row, 
          compare.phages.good.phages, 
          label_size = 12, ncol=1,
          rel_heights = c(1,1,1.5),
          labels =c('','','F'))

# Figure 4B, some only
ggsave("Figures/phages_found_metagenomes.png", width = 8.5, height = 11)

head(vibrant_results)

vibrant_results$Sample <- str_replace(vibrant_results$Sample.Name, "ME_","")
vibrant_results$Sample <- str_replace(vibrant_results$Sample, "_Bacteria.*","")
vibrant_results$Sample <- str_replace(vibrant_results$Sample, "_Virus.*","")

vibrant_results <- vibrant_results %>% gather("phageType","count",Lytic:Lysogenic)

samples.metadata <- read.csv("Data/Samples-mendota.csv")

vibrant_results <- left_join(vibrant_results, samples.metadata, by=c("Sample"="Sample.Name"))

