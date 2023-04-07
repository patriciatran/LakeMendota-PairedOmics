library(tidyverse)

checkm_combined <- read.csv("Data/checkm_output.csv", header=TRUE)
list.drep.mags <- read.csv("Data/list_dereplicated_mags", header=FALSE)



final_mags <- checkm_combined %>% filter(genome %in% list.drep.mags$V1)
colnames(final_mags)[1] <- "user_genome"

bact.gtdbtk <- read.table("Data/gtdbtk.bac120.summary.tsv", sep = "\t", header=TRUE)

arc.gtdbtk <- read.table("Data/gtdbtk.ar122.summary.tsv", sep = "\t", header=TRUE)

bact.gtdbtk <- bact.gtdbtk %>% select(user_genome, classification)

arc.gtdbtk <- arc.gtdbtk %>% select(user_genome, classification)

final_mags$user_genome <- str_replace(final_mags$user_genome, ".fasta","")


str(bact.gtdbtk$red_value)

str(final_mags)
dim(final_mags)

final_mags <- left_join(final_mags, bact.gtdbtk)

final_mags <- left_join(final_mags, arc.gtdbtk, by=c("user_genome"))

final_mags$classification <- ""

for (i in 1:nrow(final_mags)){
  if (!is.na(final_mags$classification.x[i])){
    final_mags$classification[i] <- final_mags$classification.x[i]
  }
  else{
    final_mags$classification[i] <- final_mags$classification.y[i]
    print(final_mags$classification.y[i])
  }
  print(i)
}

final_mags <- final_mags %>% select(-classification.x, -classification.y)

final_mags1 <- final_mags %>% separate(classification, into =c("Domain","Phylum", "Class","Order","Family","Genus","Species"), 
                        sep = ";", remove=FALSE)

final_mags1 %>% group_by(Domain, Phylum, Class) %>% tally() %>% arrange(desc(n)) %>%
  ggplot(aes(y=Class, x=n))+
  geom_col()+
  facet_grid(Domain+Phylum~.,scales="free", space = "free", switch="y")+
  theme_bw()+
  theme(strip.text.y.left = element_text(angle=0),
        strip.placement = "outside",
        strip.background = element_rect(fill="white"))

#ggsave("Figures/BarPlot_dRep_MAGS.pdf", width = 11, height = 8.5)

final_mags1 %>% ggplot(aes(x=contamination, y=completeness, col=Phylum))+
  geom_point(alpha=0.6)+
  theme_bw()
