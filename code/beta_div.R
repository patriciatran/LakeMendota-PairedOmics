library("phyloseq")
library("ggplot2")
library("dplyr")
library("ggpubr")
library("Matrix")
library("reshape2")
library("vegan")

mag_table <- read.table("Data/431MAGS-vs-Metagenomes-CoverM.tsv", sep="\t", header= TRUE)
sample_table <- read.csv("Data/Samples-mendota.csv")
taxa_table <- read.table("Data/MAG_taxonomy.tsv", sep="\t", header=TRUE)

mag_table <- mag_table %>% gather("READS","abund", 2:17)
mag_table$READS <- str_replace(mag_table$READS, "X", "")
mag_table$READS <- str_replace(mag_table$READS, ".fastq.*", "")


reads.match <- read.table("Data/fastq_names_match.tsv", sep="\t", header=TRUE) 

mag_table <- left_join(mag_table, reads.match)
mag_table$Sample <- str_replace(mag_table$SAMPLE.NAME, "ME_","")
mag_table$Sample <- str_replace(mag_table$Sample, "_Bacteria.*","")
mag_table$Sample <- str_replace(mag_table$Sample, "2020.","2020-")
mag_table$Sample <- str_replace(mag_table$Sample, "07.","07-")
mag_table$Sample <- str_replace(mag_table$Sample, "08.","08-")
mag_table$Sample <- str_replace(mag_table$Sample, "09.","09-")
mag_table$Sample <- str_replace(mag_table$Sample, "10.","10-")

unique(mag_table$Sample)
unique(sample_table$Sample.Name)

mag_matrix <- mag_table %>% select(Genome, Sample, abund) %>% filter(Genome != "unmapped")
mag_matrix <- mag_matrix %>% spread(Sample, abund)
mag.rows <- mag_matrix$Genome

mag_matrix <- mag_matrix[,-1]
rownames(mag_matrix) <- mag.rows

# Turn things  into phyloseq objects:
MAG <- otu_table(mag_matrix, taxa_are_rows = TRUE)
rownames(sample_table) <- sample_table$Sample.Name
sample_table[-1]

# Sample PS object:
SAMPLE <- sample_data(object = sample_table[-1])

# Taxa PS object:
taxas <- taxa_table$Genome
rownames(taxa_table) <- taxas
taxa_table <- taxa_table[,-1]
taxa_table
TAXA = tax_table(as.matrix(taxa_table))
TAXA[1:5]

# Put them all together into 1 phyloseq object:
ps.object = phyloseq(MAG, TAXA, SAMPLE)

ps.object


# Test a plot
plot_bar(ps.object, fill = "Domain")

# Relativize
relab_mag = transform_sample_counts(mag_matrix.ps, 
                                       function(x) x / sum(x) * 100) 

relab_mag

# Rename PS object:
ps.object.norm <- phyloseq(relab_mag, TAXA, SAMPLE)

# make some plots:
plot_bar(ps.object.norm, fill = "Domain")
plot_heatmap(ps.object.norm)


## BETA DIVERSITY ##
## beta diversity defined as gamma/alpha - 1:
## alpha is the average no. of species in a group, and gamma is the
## total number of species in the group

dim(SAMPLE)
dim(MAG)

# Need to transpose the table so that the length of SAMPLE and MAG are the same.
(alpha <- with(SAMPLE, tapply(specnumber(t(MAG)), Oxygen.level, mean)))
(gamma <- with(SAMPLE, specnumber(t(MAG), Oxygen.level)))
gamma/alpha - 1
