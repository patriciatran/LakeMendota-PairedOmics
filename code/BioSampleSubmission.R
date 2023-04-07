# Patricia Tran
# Script to create a TSV file to upload to NCBI for a batch genome upload
# Only 400 genomes can be uploaded at the time 

library(tidyverse)
library(data.table)

biosamples <- read.table("Data/table_download.tsv", sep="\t", header=TRUE)
sample.metadata <- read.csv("Data/Samples-mendota.csv")
MAGS <- read.table("Data/MAG_Characteristics.tsv", sep="\t", header=TRUE)
GaID <- read.csv("Data/GaID_match.csv")
 # Get the number of bases sequences in each sample.

JGI.report <- read.csv("~/Downloads/spReport (6).csv")
colnames(JGI.report)
JGI.report <- JGI.report %>% select(Sequencing.Project.Name, Total.Bases) %>% filter(Total.Bases != "NA")
JGI.report$Sequencing.Project.Name <- str_replace(JGI.report$Sequencing.Project.Name, "ME_","")
JGI.report$Sequencing.Project.Name <- str_replace(JGI.report$Sequencing.Project.Name, "Virus_.*","")

JGI.report$Sequencing.Project.Name <- str_replace(JGI.report$Sequencing.Project.Name, "_Bacteria.*","")



head(GaID)

GaID <- GaID %>% select(SAMPLE.NAME, GaID) 
GaID$SAMPLE.NAME <- str_replace(GaID$SAMPLE.NAME, "ME_","")
GaID$SAMPLE.NAME <- str_replace(GaID$SAMPLE.NAME, "_Virus.*","")
GaID$SAMPLE.NAME <- str_replace(GaID$SAMPLE.NAME, "_Bacteria.*","")
GaID

head(biosamples)
sample.metadata
biosamples$Sample.Name <- str_replace(biosamples$BioSample.name, "ME_","")

head(biosamples)
head(sample.metadata)
biosamples.with.metadata <- full_join(biosamples, sample.metadata, by="Sample.Name")
biosamples.with.metadata <- full_join(biosamples.with.metadata, GaID, by=c("Sample.Name"="SAMPLE.NAME"))


MAGS$Sample.Origin <- str_replace(MAGS$genome, "_.*","")
head(MAGS)

MAGS.with.metadata <- left_join(MAGS, biosamples.with.metadata, by=c("Sample.Origin"="GaID") )
MAGS.with.metadata$GenomeFile <- paste0(MAGS.with.metadata$genome,".fasta")
MAGS.with.metadata$AssemblyMethod <- "SPAdes"
MAGS.with.metadata$AssemblyVersion <- "v.3.13.0"



MAGS.with.metadata$SequencingTechnology <- "Illumina NovaSeq"

head(JGI.report)
head(MAGS.with.metadata)
MAGS.with.metadata <- left_join(MAGS.with.metadata, JGI.report, by=c("Sample.Name"="Sequencing.Project.Name"))

# To get the coverage, we used CoverM genome -m reads_per_bases to get the reads mapped at each based divided by the lenght of the genome. 
# But we need to get that value for the right sample.
# TMPDIR=. coverm genome --interleaved /slowdata/Reads/Paired_Bact_Vir_Mendota/Metagenomes/*.fastq -d /storage1/data10/MAGS_Final/ -m reads_per_base --min-covered-fraction 0 -o 431MAGs_coverM_reads_per_bases.tsv -x fasta -t 15


#MAGS.with.metadata$Coverage = MAGS.with.metadata$Total.Bases/MAGS.with.metadata$totalLength
#coverm <- read.table("Data/431MAGs_coverM_reads_per_bases.tsv", sep="\t", header=TRUE)

# We ran inStrain quick_profile (coverM), concatenate (rbind) all the files, and then use spread to put it in a format where each row is a MAG and each column is a sample.
coverm <- read.csv("Data/431MAGs_instrain_quickprofile.csv",
                     header=TRUE)

coverm <- coverm %>% filter(genome != "genome")
colnames(coverm)[6] <- "READS"

#View(coverm) 

#coverm <- coverm %>% gather("READS","Coverage", 2:17)
coverm$READS <- str_replace(coverm$READS, ".*_52","52")
coverm$READS <- str_replace(coverm$READS, ".fastq.*","")
head(coverm)
fastq.match <- read.table("Data/fastq_names_match.tsv", sep="\t", header=TRUE) %>% 
  select(READS, SAMPLE.NAME) %>%
  filter(grepl("Bacteria", SAMPLE.NAME))

fastq.match$SAMPLE.NAME <- str_replace(fastq.match$SAMPLE.NAME, "ME_", "")
fastq.match$SAMPLE.NAME <- str_replace(fastq.match$SAMPLE.NAME, "_Bacteria.*", "")
fastq.match$SAMPLE.NAME <- str_replace(fastq.match$SAMPLE.NAME, ".METAGENOME", "-METAGENOME")

head(fastq.match)
head(coverm)

unique(coverm$READS)
coverm$READS <- str_replace(coverm$READS, "-",".")
coverm$READS <- str_replace(coverm$READS, "-",".")
unique(fastq.match$READS)

coverm <- left_join(coverm, fastq.match)
head(coverm)

coverm <- coverm[,c(1,4,6,7)]
colnames(coverm)[4] <- "Sample.Name"
colnames(coverm)[1] <- "genome"
head(coverm)
coverm$Sample.Name <- str_replace(coverm$Sample.Name, "\\.","-")
coverm$Sample.Name <- str_replace(coverm$Sample.Name, "\\.","-")
coverm$genome <- str_replace(coverm$genome, ".fasta","")

MAGS.with.metadata$Sample.Name
coverm$Sample.Name

head(coverm)

head(MAGS.with.metadata)
test1 <- left_join(MAGS.with.metadata, coverm) %>% distinct()

head(MAGS.with.metadata$genome)
head(coverm$genome)

unique(MAGS.with.metadata$Sample.Name)
unique(coverm$Sample.Name)

test1$coverage
#colnames(test1)
#test1$READS

test1$AssemblyName <- paste0("UWMLM-",
                                    toupper(substr(test1$Phylum, 4,6)),
                                    "-",
                                    toupper(substr(test1$Class, 4,6)),
                             "-",
                             toupper(substr(test1$Order, 4,6)),
                             "-",
                             toupper(substr(test1$Genus, 4,6)),
                             "-",
                             test1$genome,
                                    "-v1.0")

#### CREATE BIOSAMPLE TAXA ####
BioSamples.Submission <- test1 %>% select(genome, classification, Date.YYYY.MM.DD, Sample.Depth.meters) %>%
  mutate(BioProject_accession = "PRJNA758276",
         isolate = "freshwater pelagic",
         ENVO = 'ENVO:00000021', # freshwater lake,
    Local_ENVO = 'ENVO:06105011',
    Env_medium_ENVO = 'ENVO:04000007',
    Geography = 'USA:Wisconsin',
    Lat_Long = '43.1097 N 89.4206 W',
    PhysicalSource = paste0(genome, ' binned from freshwater Lake Lake Mendota WI water column')
  ) %>% 
  distinct()

write.table(BioSamples.Submission, "Data/BioSamples.tsv",
            sep="\t",
            quote = FALSE, 
            row.names = FALSE)    

## Submit Genomes:
MAG.biosamples <- read.table("Data/BioSampleObjects.txt", sep="\t", header=TRUE)[,1:2]
# Only get the columns needed and sort the names alphabetically.

submission.table <- test1 %>% select(AssemblyName, genome, 
                                     AssemblyMethod, AssemblyVersion, 
                                     SequencingTechnology, coverage,
                                     GenomeFile
) %>% left_join(MAG.biosamples, by=c("genome"="Sample.Name")) %>%
  arrange(genome) 
  


test2 <- submission.table %>% distinct()
test2$Coverage <- paste0(substr(test2$coverage, 1,4),"x")

test2 <- test2 %>% arrange(genome)

#test2
length(unique(test2$AssemblyName))

write.table(test2[1:300,], "Data/Submission_Genomes_PartI_newCov.tsv",
            sep="\t",
            row.names = FALSE,
            quote = FALSE)


write.table(test2[301:431,], "Data/Submission_Genomes_PartII_newCov.tsv",
            sep="\t",
            row.names = FALSE,
            quote = FALSE)

write.table(test2,
            "Data/Submission_newCoverage.tsv",
            sep="\t",
            row.names = FALSE,
            quote = FALSE)

## merge physical metadata with MAG data:
BioSamples16_accesions <-read.table("Data/BioSampleAccession16PhysicalSamples.tsv", sep="\t", header=TRUE)
head(BioSamples16_accesions) 
BioSamples16_accesions <- BioSamples16_accesions %>% select(accession, sample_name)
BioSamples16_accesions$sample_name <- str_replace(BioSamples16_accesions$sample_name, "ME_","")
colnames(BioSamples16_accesions) <- c("phys_sample_accession","phys_sample_name")


BioSamplesMAG_accession <- read.table("Data/BioSamples431MAGsAccession.tsv", header=TRUE,
                                      sep="\t") %>% select(accession, sample_name, bioproject_accession)

BioSamplesMAG_accession$GaID <- str_replace(BioSamplesMAG_accession$sample_name, "_.*","")

BioSamplesMAG_accession <- BioSamplesMAG_accession %>% left_join(GaID)

BioSamplesMAG_accession_with_phys_sample <- BioSamplesMAG_accession %>% left_join(BioSamples16_accesions, 
                                                                                  by=c("SAMPLE.NAME" = "phys_sample_name"))


colnames(BioSamplesMAG_accession_with_phys_sample) <- c("BioSample MAG accession",
                                                        "sample_name_MAG",
                                                        "bioproject_accession","JGI GaID",
                                                        "physical sample name",
                                                        "physical sample BioSample Accession")


write.table(BioSamplesMAG_accession_with_phys_sample,
            "Data/Data_to_email_NCBI_MAG_linked_to_phys_samples.tsv",
            sep="\t",
            row.names = FALSE,
            quote = FALSE)
