phages.rpkm.expr <- read.table("Data/phages_rpkm_normalized_matrix.tsv", sep="\t", header=TRUE)

total_phages <- nrow(phages.rpkm.expr)

samples <- colnames(phages.rpkm.expr)[2:17]

index=1

checkv <- read.table("Data/Phages_all_qualities_and_sizes_with_taxonomy.tsv", 
                     sep="\t",
                     header=TRUE)

phage.type <- checkv %>% select(contig_id, type_by_vibrant)


colnames(phages.rpkm.expr)
colnames(phage.type)

phages.rpkm.expr <- left_join(phages.rpkm.expr, phage.type, 
                              by=c("Genome"="contig_id"))


percent_expressed <- ""

for (i in 2:ncol(phages.rpkm.expr)){
  phages.sample <- phages.rpkm.expr[,i]
  expressed_once <- length(which(phages.sample > 0))
  format(expressed_once/total_phages*100, digits=3)
  #colnames(phages.rpkm.expr)[i]
  percent_expressed[index] <- format(expressed_once/total_phages*100, digits=3)
  percent_expressed

  index = index + 1
  
}



samples
percent_expressed

phages.percentage <- as.data.frame(matrix(c(samples, 
                                            as.numeric(percent_expressed)
                                            ), ncol=2, byrow = FALSE))

phages.percentage$V1 <- str_replace(phages.percentage$V1 ,"X","")
phages.percentage$V1 <- str_replace(phages.percentage$V1 ,"\\.","-")
phages.percentage$V1 <- str_replace(phages.percentage$V1 ,"\\.","-")

phages.percentage$V1 

colnames(phages.percentage) <- c("Sample.Name","Percentage")
#sample.metadata
library(tidyverse)
library(lubridate)

sample.metadata <- read.csv("Data/Samples-mendota.csv")

phages.percentage <- left_join(phages.percentage, sample.metadata)

cols_oxygen <- c("oxic"="#48AA72",
                 "oxycline"="#red",
                 "anoxic"="#3E1E53")

max(phages.percentage$Percentage)
min(phages.percentage$Percentage)

ggplot(phages.percentage, aes(x=as.Date(Date.YYYY.MM.DD), 
                              y=as.numeric(Percentage), 
                              col=Oxygen.level,
                              shape=as.factor(Sample.Depth.meters)
))+
  geom_point(size=4)+
  #geom_line()+
  scale_color_manual(values = cols_oxygen)+
  theme_bw()+
  ylab("Percentage of high quality phages\n that are expressed")+
  xlab("Sample date")+
  scale_x_date(breaks = unique(ymd(sample.metadata$Date.YYYY.MM.DD)),
               date_labels = "%b %d")+
  theme_bw()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle=0),
        text = element_text(size=12))

# Figure4D
ggsave("Figures/Percent_of_active_phages_over_time.pdf", width = 11, height = 6)



phages.rpkm.expr %>%  gather("sample","expression", 2:17) %>% 
  mutate(Activity = ifelse(expression == 0, yes="Inactive",no="Active")) %>%
  group_by(sample, Activity, type_by_vibrant) %>% tally() %>%
  ggplot(aes(x=sample, y=n, fill=type_by_vibrant))+
  geom_col()+
  facet_grid(Activity~., scales="free")+
  theme_bw()+
  theme(strip.background = element_blank(),
        axis.text.x = element_text(angle=90, hjust=1),
        panel.grid.minor = element_blank(),
        strip.text.y = element_text(angle=0),
        text = element_text(size=12),
        panel.grid = element_blank())

ggsave("Figures/Phages_active_inactive_type.pdf", width = 10, height = 7)
  

## Assign links:
links <- read.csv("Data/Host_prediction_to_genome_m90.csv")

phage.rpkpm.express.host <- left_join(phages.rpkm.expr, 
                                      links, 
                                      by=c("Phage" = "Virus"))

write.table(phage.rpkpm.express.host, "Data/phage.rpkm.expression.with.host.tsv",
            sep='\t', row.names = FALSE)

# Get the rank for each phage across each sample
dim(phages.rpkm.expr)
phages.rank <- matrix(nrow=nrow(phages.rpkm.expr), 
                      ncol=ncol(phages.rpkm.expr))
colnames(phages.rank) <- colnames(phages.rpkm.expr)
phages.rank[,1] <- phages.rpkm.expr[,1]
phages.rank <- as.data.frame(phages.rank)


for (i in 2:ncol(phages.rpkm.expr)){
  sample.name.to.get <- colnames(phages.rpkm.expr)[i]
  
  sample.name.to.get
  test1 <- phages.rpkm.expr[c(1,i)]
  test1 <- test1[order(test1$X2020.07.24_15m,decreasing=TRUE),]
  
  rank=1
  for (phage in 1:nrow(test1)){
    phage.name <- test1$Phage[phage]
    phages.rank[phages.rank[phage.name],]
    
  }
  

}





