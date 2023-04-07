## Patricia Tran
## Feb 27 2023
## Description: This script consolidates the genome annotation data
## and gene-level expression data to Create Figure 7. This also 
## incorporates taxonomic information, and whether the MAG was 
## phage-associated or not.

metabolic.table <- read.table("Data/METABOLIC_table.tsv", 
                              sep="\t", 
                              header=TRUE,
                              quote = "")

library(tidyverse)

metabolic.table <- metabolic.table %>% gather("MAG","Presence", 11:1303)

# Remove the rows with the word "numbers"
metabolic.table1 <- metabolic.table %>% filter(grepl(pattern = "presence", x = MAG))
metabolic.table1$MAG <- str_replace(metabolic.table1$MAG, ".Hmm.presence","")
taxonomy <- read.table("Data/MAG_taxonomy.tsv", sep="\t", header=TRUE)         

metabolic.table1 <- left_join(metabolic.table1, taxonomy, by=c("MAG"="Genome"))


#### RENAME TAXA ####

metabolic.table1 <- metabolic.table1 %>% mutate(Taxa_to_plot = ifelse(Phylum == "p__Proteobacteria",
                                                                      yes=Class,
                                                                      no=Phylum))


metabolic.table1 <- metabolic.table1 %>% mutate(Taxa_to_plot2 = ifelse(Taxa_to_plot == "c__Gammaproteobacteria",yes=Order,
                                                                       no=Taxa_to_plot))



# My color scheme:
cols2 <- c("p__Verrucomicrobiota"="#AF1F80",
           "p__Sumerlaeota"="grey",
           "p__Planctomycetota"="#FF0000",
           "p__Patescibacteria"="grey",
           "p__Nanoarchaeota"="grey",
           "p__Myxococcota"="grey",
           "p__Krumholzibacteriota"="grey",
           "p__Hydrogenedentota"="grey",
           "p__Gemmatimonadota"="grey",
           "p__Firmicutes"="#CFE2F3",
           "p__Firmicutes_A"="#CFE2F3",
           "p__Firestonebacteria"="grey",
           "p__Fibrobacterota"="grey",
           "p__Desulfobacterota"="#B4A7D6",
           "p__Desulfobacterota_F"="#B4A7D6",
           "p__Dependentiae"="grey",
           "p__Cyanobacteria"="#8BAF7B",
           "p__Chloroflexota"="#6100FF",
           "p__Chlamydiota"="grey",
           "p__Campylobacterota"="grey",
           "p__Bdellovibrionota"="grey",
           "p__Bdellovibrionota_C"="grey",
           "p__Bacteroidota"="grey",
           "p__Armatimonadota"="grey",
           "p__Actinobacteriota"="#DEA955",
           "p__Acidobacteriota"="grey",
           "o__Thiomicrospirales"="grey",
           "o__Steroidobacterales"="grey",
           "o__Pseudomonadales"="grey",
           "o__Methylococcales"="#2F83D4",
           "o__GCA-2729495"="#2FB5D4",
           "o__Burkholderiales"="#2BCFF4",
           "c__Alphaproteobacteria"="#EA9999")


MAGs.with.phages <-read.csv("Data/MAGs-with-a-phage-manual-or-iPHOP.csv")
unique(MAGs.with.phages$HostGenome)

MAGS.with.phages1 <-MAGs.with.phages %>% select(HostGenome) %>% distinct() %>% 
  mutate(HasPhage = "Yes")


head(metabolic.table1)


# Save file:
write.table(metabolic.table1, "Data/METABOLIC_Table.clean.tsv", sep="\t", quote = FALSE, row.names = FALSE)


unique(metabolic.table1$Category) # get nitrogen cycling

#### DIAGRAMS FOR EACH CYCLE####
## Nitrogen cycle:
nitrogen.cycling <- metabolic.table1 %>% filter(Category == "Nitrogen cycling") 
nitrogen_cycling_phage <- left_join(nitrogen.cycling, MAGS.with.phages1, by=c("MAG"="HostGenome"))


nitrogen_cycling_phage$HasPhage <- replace_na(nitrogen_cycling_phage$HasPhage, "No")

nitrogen.plot <- nitrogen_cycling_phage %>% filter(Presence == "Present") %>%
  group_by(Function, Gene.abbreviation,Presence, HasPhage) %>% tally() %>%
  ggplot(aes(x=n, y=Gene.abbreviation, fill=HasPhage))+
  geom_bar(stat ="identity", colour="black")+
  scale_fill_manual(values=c("white","grey"))+
  facet_grid(Function~., scales="free_y")+
  theme_bw()+
  theme(strip.text.y = element_text(angle=0),
        legend.position = "bottom",
        strip.background = element_blank(),
        panel.grid = element_blank())+
  xlab("Total count of MAGs with this gene")+
  geom_text(aes(label=n),position="stack",hjust=1, vjust=0.5)
nitrogen.plot
## Sulfur
sulf.cycling <- metabolic.table1 %>% filter(Category == "Sulfur cycling") 
sulfur_cycling_phage <- left_join(sulf.cycling, MAGS.with.phages1, by=c("MAG"="HostGenome"))


sulfur_cycling_phage$HasPhage <- replace_na(sulfur_cycling_phage$HasPhage, "No")

sulfur.plot <- sulfur_cycling_phage %>% filter(Presence == "Present") %>%
  group_by(Function, Gene.abbreviation,Presence, HasPhage) %>% tally() %>%
  ggplot(aes(x=n, y=Gene.abbreviation, fill=HasPhage))+
  geom_bar(stat ="identity", colour="black")+
  scale_fill_manual(values=c("white","grey"))+
  facet_grid(Function~., scales="free_y")+
  theme_bw()+
  theme(strip.text.y = element_text(angle=0),
        legend.position = "bottom",
        strip.background = element_blank(),
        panel.grid = element_blank())+
  xlab("Total count of MAGs with this gene")+
  geom_text(aes(label=n),position="stack",hjust=1, vjust=0.5)

sulfur.plot  

sulfur_cycling_phage %>% filter(Function == "Sulfate reduction" &
                                  Presence == "Present" &
                                  HasPhage == "Yes") %>% 
  select(MAG, Phylum, Gene.abbreviation) 

library(ggpubr)
ggarrange(nitrogen.plot, sulfur.plot, nrow=2,
          labels=c("A","B"), common.legend = TRUE)

ggsave("Figures/Nitrogen_Sulfur_phage_vs_nophages.pdf", width = 10, height = 7)

#### METHANE, N and S for PHAGE-ASSOCIATED MAGS ####
table1 <- metabolic.table1 %>% filter(Category %in% c("Methane metabolism", "Nitrogen cycling", "Sulfur cycling")) %>% 
  left_join(MAGS.with.phages1, by=c("MAG"="HostGenome")) %>%
  filter(Presence == "Present") %>%
  group_by(Category, Function, Gene.abbreviation,Presence, HasPhage) %>% tally() 

table1$HasPhage <- replace_na(table1$HasPhage, "No")
  
ggplot(table1, aes(x=n, y=Gene.abbreviation, fill=HasPhage), width=0.5)+
  geom_bar(stat ="identity", colour="black")+
  scale_fill_manual(values=c("white","grey"))+
  facet_grid(Category+Function~., scales="free", space = "free_y")+
  theme_bw()+
  theme(strip.text.y = element_text(angle=0, vjust=1),
        legend.position = "bottom",
        strip.background = element_blank(),
        panel.grid = element_blank())+
  xlab("Total count of MAGs with this gene")+
  geom_text(aes(label=n),position="stack",hjust=1, vjust=0.5, size=3)

ggsave("Figures/Methane_S_N.pdf", width = 8.5, height = 11)


####  GENE EXPRESSION OF THESE CYCLES OF INTEREST ####
library(data.table)

# Load all the gene RPKM
gene.coverage <- fread("Data/All_gene_collections_transcript_coverage_cat.txt")
gene.headers <- fread("Data/all_hmm_collection.headers-1.txt", header = FALSE)[,-1]

head(gene.coverage)
head(gene.headers)


colnames(gene.coverage) <- c("contigName","RPKM","Sample")
colnames(gene.headers) <- c("HMMFileName","contigName")
# Get only those that are expressed:
gene.coverage.expressed <- gene.coverage %>% filter(RPKM > 0)
dim(gene.coverage.expressed)

# There are 3,308,642 genes expressed across all the MAGs

# Now we  don't care when the sample was expressed, only if it WAS expressed, so let's remved the sample info 
genes.expressed.once <- gene.coverage.expressed %>% select(-Sample) %>% distinct()

head(gene.headers)
head(genes.expressed.once)


gene.expressed.once.with.names <- left_join(genes.expressed.once, gene.headers, by="contigName")

#head(gene.coverage)
#head(genes.expressed.once)
head(gene.expressed.once.with.names)

# Only keep those that have an HMM found:
gene.expressed.once.with.names <- gene.expressed.once.with.names %>%  filter(!is.na(HMMFileName))

# Split the contigName by MAG:
gene.expressed.once.with.names <- gene.expressed.once.with.names %>%  separate(contigName, into=c("MAG","Contig"), 
                                            sep="~~", remove = FALSE)

## Create a table that matches the .HMM file names to the genes:
head(metabolic.table1)

head(metabolic.table1[1:5,1:5])

# We want to split thhe column Hmmm file by the comma and create a row for it:
metabolic.table1_sep_hmm <- separate_rows(data = metabolic.table1, sep=", ")

nrow(metabolic.table1_sep_hmm)
head(metabolic.table1_sep_hmm )
unique(metabolic.table1_sep_hmm$Hmm.file)
head(gene.expressed.once.with.names)


## Now for each row there is a row, get the Hmmname and the MAG name, and then
# if there is a corresponding row in gene.expressed.once.with.names, then it means the gene was active.



pairs <- gene.expressed.once.with.names %>% select(MAG, HMMFileName) %>% distinct()


fewer_pairs <-pairs %>% filter(HMMFileName %in% unique(metabolic.table1_sep_hmm$Hmm.file))
fewer_pairs$Expressed <- "Active"

metabolic.table1_sep_hmm_with_expression <-left_join(metabolic.table1_sep_hmm, fewer_pairs, 
          by=c("MAG"="MAG",
               "Hmm.file" ="HMMFileName"))


metabolic.table1_sep_hmm_with_expression$Expressed <- replace_na(metabolic.table1_sep_hmm_with_expression$Expressed, "Inactive")

summary(metabolic.table1_sep_hmm_with_expression)

unique(metabolic.table1_sep_hmm_with_expression$Expressed)

write.table(metabolic.table1_sep_hmm_with_expression, "Data/SuppTable6-ActiveInactiveGenes.tsv",
            sep="\t",
            quote=FALSE,
            row.names = FALSE)

#### READY TO PLOT WITH THE PHAGES ####

data_to_plot <-metabolic.table1_sep_hmm_with_expression %>% 
  #filter(Function %in% c("Wood Ljungdahl pathway")) %>%
  filter(Category %in% c("Methane metabolism", "Nitrogen cycling", "Sulfur cycling")) %>% 
  left_join(MAGS.with.phages1, by=c("MAG"="HostGenome")) %>%
  filter(Presence == "Present" & HasPhage == "Yes") %>%
  group_by(Category, Function, Gene.abbreviation,
           Taxa_to_plot2, MAG, HasPhage, Expressed) %>% tally()

colnames(data_to_plot)
colnames(gene.expressed.once.with.names)

length(unique(data_to_plot$MAG))


#### PLOT IT ####
ggplot(data_to_plot, aes(y=MAG, x=Gene.abbreviation, 
                         shape=Expressed,
                         col=Category
                         ))+
  geom_point(size=3)+
  facet_grid(Taxa_to_plot2~Category+Function, scales="free",space="free")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_text(angle=0),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        text = element_text(color="black"),
        legend.position="top")+
  scale_shape_manual(values=c(16,1))+
  scale_color_manual(values=c("#5296a5","#7dce82",
                              "#961d4e"))



ggsave("Figures/methane_sulfur_nitrogen_MAGS_with_phage_with_expression_pattern.pdf", 
       width = 11, height = 8)

## Plot all:
metabolic.table1_sep_hmm_with_expression %>% 
  #filter(Function == "Wood Ljungdahl pathway") %>%
  filter(Category %in% c("Methane metabolism", "Nitrogen cycling", "Sulfur cycling")) %>% 
  left_join(MAGS.with.phages1, by=c("MAG"="HostGenome")) %>%
  filter(Presence == "Present") %>%
  group_by(Category, Function, Gene.abbreviation,
           Taxa_to_plot2, MAG, HasPhage, Expressed) %>% tally() %>%
  ggplot(aes(y=MAG, x=Gene.abbreviation, 
                           shape=Expressed,
                           col=Category
  ))+
  geom_point(size=3)+
  facet_grid(HasPhage+Taxa_to_plot2~Category+Function, scales="free",space="free")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_text(angle=0),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        text = element_text(color="black"),
        legend.position="top")+
  scale_shape_manual(values=c(16,1))+
  scale_color_manual(values=c("#5296a5","#7dce82",
                              "#961d4e"))

#ggsave("Figures/NSM.activity.MAGs.pdf", width = 10, height = 10)
  


#### OTHER ####
unique(MAGs.with.phages$HostGenome)

phages.per.mag <- MAGs.with.phages %>% select(HostGenome, X) %>% distinct() %>% 
  group_by(HostGenome) %>% tally()
colnames(phages.per.mag) <- c("MAG","Phages per MAG")


test2 <- metabolic.table1 %>% filter(Category %in% c("Methane metabolism", "Nitrogen cycling", "Sulfur cycling")) %>% 
  left_join(MAGS.with.phages1, by=c("MAG"="HostGenome")) %>%
  left_join(phages.per.mag) %>%
  filter(Presence == "Present" & HasPhage == "Yes") %>%
  select(MAG, Taxa_to_plot2, `Phages per MAG`) %>% distinct() %>% arrange(Taxa_to_plot2)

ggplot(test2, aes(x=`Phages per MAG`, y=MAG))+
  geom_col()+
  facet_grid(Taxa_to_plot2~., scales = "free", space = "free")+
  theme_bw()+
  theme(strip.background = element_blank(),
        strip.text.y = element_text(angle=0),
        panel.grid = element_blank(),
        axis.text.x = element_text(angle=45, hjust=1),
        text = element_text(color="black"),
        legend.position="none")

## Table format:
metabolic_table_phage.host <- metabolic.table1 %>% filter(Category %in% c("C1 metabolism", "Nitrogen cycling", "Sulfur cycling")) %>% 
  left_join(MAGS.with.phages1, by=c("MAG"="HostGenome")) %>%
  filter(Presence == "Present") %>%
  group_by(Function, Gene.abbreviation,Presence, HasPhage) %>% tally() %>%
  spread(HasPhage, n)

write.table(metabolic_table_phage.host, "Data/Metabolic_Table_Phage_Host_N_and_S.tsv", sep="\t",
            row.names = FALSE, quote = FALSE)


NScycle_by_MAG_with_phage <- metabolic.table1 %>% filter(Category %in% c("Nitrogen cycling", "Sulfur cycling")) %>% 
  left_join(MAGS.with.phages1, by=c("MAG"="HostGenome")) %>%
  filter(Presence == "Present" & HasPhage == "Yes") %>%
  group_by(Function, Gene.abbreviation, Taxa_to_plot2, MAG) %>% tally() %>%
  spread(MAG, n)

write.table(NScycle_by_MAG_with_phage, "Data/Metabolic_Table_Phage_Host_N_and_S_by_MAG.tsv", sep="\t",
            row.names = FALSE, quote = FALSE)

## Anther plot:

nitrogen_cycling_phage %>% filter(Presence == "Present") %>%
  group_by(Function, Gene.abbreviation,Presence, Taxa_to_plot2, HasPhage) %>% tally() %>%
  ggplot(aes(x=n, y=Gene.abbreviation, fill=Taxa_to_plot2))+
  geom_bar(stat ="identity", colour="black")+
  #scale_fill_manual(values=c("white","grey"))+
  scale_fill_manual(values = cols2)+
  facet_grid(Function~., scales="free_y")+
  theme_bw()+
  theme(strip.text.y = element_text(angle=0),
        legend.position = "bottom",
        strip.background = element_blank(),
        panel.grid = element_blank())+
  xlab("Total count of MAGs with this gene")

sulfur_cycling_phage %>% filter(Presence == "Present") %>%
  group_by(Function, Gene.abbreviation,Presence, Taxa_to_plot2, HasPhage) %>% tally() %>%
  ggplot(aes(x=n, y=Gene.abbreviation, fill=Taxa_to_plot2))+
  geom_bar(stat ="identity", colour="black")+
  #scale_fill_manual(values=c("white","grey"))+
  scale_fill_manual(values = cols2)+
  facet_grid(Function~., scales="free_y")+
  theme_bw()+
  theme(strip.text.y = element_text(angle=0),
        legend.position = "bottom",
        strip.background = element_blank(),
        panel.grid = element_blank())+
  xlab("Total count of MAGs with this gene")
