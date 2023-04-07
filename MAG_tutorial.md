## The data is released to the public and is under 2 year embargo.

| JGI Proposal | 506328 |
| --- | --- |
| Number of bacterial samples | 16 |
| Number of viral samples | 14 |
| Time span  | July 2020 to October 2020 |
| Location | Deep Hole Lake Mendota |
| NCBI Project | https://dataview.ncbi.nlm.nih.gov/object/PRJNA758276 |

## 🧬 Total Metagenomes & Viromes data link on the JGI:

[https://genome.jgi.doe.gov/portal/Paivirtimeseries/Paivirtimeseries.info.html](https://genome.jgi.doe.gov/portal/Paivirtimeseries/Paivirtimeseries.info.html)

## 🧬 Data Availability

The MAGs can be accessed on BioProject PRJNA758276. Under this Bioproject can be accessed all the 431 prokaryotic metagenome-assembled genomes, and the processed 16 RNA-seq samples. The phage genomes can be accessed at 10.6084/m9.figshare.22213846. (note to reviewers: private link is [https://figshare.com/s/9502a09fdba078edda07](https://figshare.com/s/9502a09fdba078edda07) ). 

# General useful commands

Skip to the [tutorial for processing on samples](https://www.notion.so/Lake-Mendota-Paper-4ca5288191e64edfab21674f31775175)

`du --summarize --human-readable *` To show the folder sizes for each subfolder

`df -h /slowdata/` to check the size of the /slowdata folder shared by all lab members

`grep '>' -c [filename.fasta]` to count how many fasta headers a file contains

`ls -lh [folder]|wc -l` to count how many elements a folder contains

`du -sh /slowdata/databases/` to check the size of a specific folder

**[Conda](https://docs.conda.io/en/latest/)** environments refresher:

```
ENVNAME=metabat2_run && \\
conda create -n $ENVNAME python=3.7 &&\\
conda install -n $ENVNAME -c biocconda metabat2 &&\\
conda activate $ENVNAME

conda env export > dbcan_condaenv_macox.yml

conda deactivate

conda env create -f dbcan_condaenv_macox.yml

```

`nohup bash [.sh file] > [logfile] &` to run a process in the background

# How to download data

# 1- Accessing Data (applicable to UW-Madison users)

Follow the instructions on the [UW DoIT page](https://kb.wisc.edu/researchdata/108855#access) to set an account:
Go on your JGI Download page for the [whole proposal](https://genome.jgi.doe.gov/portal/pages/dynamicOrganismDownload.jsf?organism=Paivirtimeseries), click the checkbox for the folders you want to download.

Follow the instructions on the DoIT page to download the files to the server directly.

Log into the lab's server, and the JGI folders are located in `/mnt/researchdrive/LakeMendotaPatricia_PairedVirBact2020/`.

# 2- Using Globus to download the necessary files from the JGI (applicable to anyone)

The JGI already performs assemblies of the file, so it is not necessary to redo the assembly. However, I suggest to check the report tables they generate (.PDF), and look at the scaffolds using a program like [Bandage](https://rrwick.github.io/Bandage/) for added peace of mind.
If you are working with the [Bacteria](notion://www.notion.so/Lake-Mendota-Paper-4ca5288191e64edfab21674f31775175#link), please see the section below. Otherwise, skip to the [Viral](notion://www.notion.so/Lake-Mendota-Paper-4ca5288191e64edfab21674f31775175#Viral-metagenomes) section.

Login to JGI Portal > go to Favorites

1. Select all projects > by Type > Assembly > `.*contigs.fna` > download using globus
2. Select all projects > by Type > Assembly > `.*prodigal_proteins.faa` > download using globus
3. Select all projects > by Type > Assembly > `.*prodigal_genes.fna` > Download using Globus
4. Select all projects > by Type > Sequence > `.*METAGENOME.fastq.gz` > download using globus [this is likely going to be the largest file download]

It’s important to get the right file with the Ga0000 ID in the headers, and not the one with scaffold_c1_1 in the headers, because it will make it easier to match the IMG annotations later on. 

Plus, there will be times when we will need to concatenate all the scaffolds together into 1 file, and we do not want 2 headers from 2 separate samples named the same way.

All in all the files that I really needed for this project were about 60GB not including the fastq.gz reads. Those are approximately 50GB per file (unzipped)

**The total space you will need is:**
(50GB*30)+60GB ~ 1.56TB

The metatranscriptomes is a dataset on its not (not on the JGI). The metatranscriptomes (forward and reverse reads in fastq format) can be downloaded from the SRA under NCBI Project ID: **PRJNA758276**

## 3- Using the API to download files directly to the server without using Globus

This is a method that can be useful if you don't have Globus enabled on your server, for example with the GLBRC server.

Steps:

1. Follow the same steps to obtain the Download link by email as above.
2. Log in to GLBRC:
3. Set your cookies to store the JGI information:

```
curl '<https://signon.jgi.doe.gov/signon/create>' --data-urlencode 'login=USER_NAME' --data-urlencode 'password=USER_PASSWORD' -c cookies

```

Replace  USER_NAME with your username (e.g. [ptran5@wisc.edu](mailto:ptran5@wisc.edu)) with the full email extension. Replace USER_PASSWORD with your own user password.

the `-c` here means write the cookies to the file named `cookies`

1. There should now be a cookies file (`ls cookies`) in your working directory that looks like this:

```
# Netscape HTTP Cookie File
# <https://curl.se/docs/http-cookies.html>
# This file was generated by libcurl! Edit at your own risk.
#HttpOnly_signon.jgi.doe.gov    FALSE   /   FALSE   0   _caliban_session    bmg2YWVFTWRmK2ljYWZ0VmRiZGVzYUxCOGIxdW9iQ1d6ZU5GS0NRYk1NcDhmUHFLZE1QV1hORVUvSjgrRkdUcDZxc2lGUXZOTXFkSjRvRS9TQ29PMjh1MHk5UGd5OUd2U0JoQ1ZZWUsydjNaRWVvWHY0Yy9FZ2hpTGE2cDJlWDlSVG5RbStGWS9RcmVMQ0x3U2R2eDduYTBmU1ZwMmFrbnNKQkIwK1B1UDhwcnhjSXFkVXdvY0d2VnlOLzB0eW5VL1pGdnVkTklYaFBQTUFhNUpyQjVIZz09LS0zUmFheGg3VTFOOS9PUUY0ZS8wb3B3PT0%3D--121afbb54f9fb112489c494ee2c38db07dde4305
.jgi.doe.gov    TRUE    /   FALSE   0   jgi_session %2Fapi%2Fsessions%2F363a77eca058a1990c1658e770e5a055

```

1. Create an executable scripts using `nano Download_JGI_data.sh`
that looks like this:

```
cd /home/GLBRCORG/ptran5/home/GLBRCORG/ptran5/

curl '<https://genome.jgi.doe.gov/portal/ext-api/downloads/bulk/75941-51/zip>' -b cookies > FAA_assembled.zip

```

Replace the link in the single quotation marks with the Download link. to find the Download link, click your Email, and right-click "Copy link address". Whatever is next to ">" will be the name of your zip file.

Save file.

Note the `-b` option which means "use the cookies in the `cookies` file we created in Step 3.

1. Make it executable
`chmod u+x Download_JGI_data.sh`
2. Create your `condor_submit` script that looks something like this:

```
[TBD]
```

1. Submit your job:
`condor_submit Download_JGI_data.sub`

## JGI File descriptions:

- contigs.fna : the assembled scaffolds by the JGI with headers Gaxxx_xxxx that matches those on the portal
- prodigal_proteins.faa : the prodigal called proteins on the contigs.fna with the headers being GaXXX_XXX_XX (amino acids)
- prodigal_proteins.faa : the prodigal called genes on the contigs.fna with the headers being GaXXX_XXX_XX (nucleotides)
- METAGENOME.fastq.gz : the quality filter FASTQ reads used to create the assemblies. They are interleaved (only 1 fastq but with forward (1) and reverse (2) in the same file) (JGI pipeline V5)

# Working directories

We will move the Assembly files (fasta) and Reads files (fastq) to `slowdata/Reads/LakeMendotaPatricia/` as we work with them. Then after finishing to work we them, please zip them using `pigz [filename]` to save some space! Our analyses will be in the folder named `/slowdata/data2/Lake_Mendota_Patricia`.

The folder structure is:

```
/slowdata/Reads/Paired_Bact_Vir_Mendota/
|___ Metagenomes
    |___ METAGENOMES.fastq.gz
    |___ assembled_contigs/
        |___ *assembly.contigs.fasta
    |___ prodigal_res/
        |___*prodigal_proteins.ffa
        |___*prodigal_genes.fna
        |___reformatting/
            |___*prodigal_genes.reformat.fna
|___ Viromes
    |___ VIRAL-METAGENOMES.fastq.gz
    |___ assembled_contigs/
        |___ .*assembly.contigs.fasta
    |___ prodigal_res/
        |___*prodigal_proteins.ffa
        |___*prodigal_genes.fna
        |___reformatting/
            |___*prodigal_genes.fna
|__ Metatranscriptomes
   |___ METATRANSCRIPTOMES.fastq.gz

```

The analyses are going to be in `/storage1/data10/`

# Bioinformatics programs needed

Here is a list of programs that we will be using. I suggest creating a conda environment for each software.

`miniconda` is a package/software and dependency manager. I like to think of it as we enter a space, and open a door (the `conda env` / `conda activate [environment name]`) to a room that contains all the dependencies for a given program to run. Each "room" is independent of each other. We can then exit (`conda deactivate`) the room and move on. All the programs mentioned in the paper can be downloaded via `conda environments` . However, there were a couple that needed more setting up so I included extra instructions below.

- MetaWRAP:  https://github.com/bxlab/metaWRAP
    
    The way I was able to install metawrap was slightly different than the one on the website.
    
    1. Login to server. outside of your conda environment, go to your /home/username directory and pull the git directory: git clone [https://github.com/bxlab/metaWRAP.git](https://github.com/bxlab/metaWRAP.git)
    2. Carefully configure the `/home/username/metaWRAP/bin/config-metawrap` file to it points to your desired database locations (you can modify this later). Follow the database configuration guide for details. The databases are already on the server and I've attached what my config-metawrap file here.
    3. copy over the contents of `yourpath/metaWRAP/bin/` into a location already in your $PATH (for us it's going to be: /home/username/miniconda3/bin/)
    
    ```bash
    conda create -y -n metawrap-env python=2.7
    conda activate metawrap-env
    ```
    
    ```bash
    5. Now it should show metawrap-env in front of your prompt. type
    
    ```
    
    ```bash
    conda config --add channels defaults
    conda config --add channels conda-forge
    conda config --add channels bioconda
    conda config --add channels ursky
    conda install --only-deps -c ursky metawrap-mg
    
    ```
    
    ```
    6. Now metawrap should work WITHIN the environment!
    
    ```
    
    `metawrap -h` to see the manual.
    
- checkM v.1.1.3 :  [https://anaconda.org/bioconda/checkm-genome](https://anaconda.org/bioconda/checkm-genome)
    
    CheckM installation is kind of default, but you will have to download the CheckM database, and use checkm setRoot to tell it where to find it.
    
    ```
    conda create -n checkm_1.1.3 python=3.7
    conda install -n checkm_1.1.3 -c bioconda checkm=1.1.3
    conda activate checkm_1.1.3
    checkm setRoot /path/to/checkmDB
    
    ```
    
- GTDB-tk v.1.7.0
[https://anaconda.org/bioconda/gtdbtk](https://anaconda.org/bioconda/gtdbtk)
    
    ( I used v1.7.0 with database release 202)
    The best way I found to install the program was definitely in a conda environment. Firt, decide which release you want (release95, 202), the [table is here](https://ecogenomics.github.io/GTDBTk/installing/index.html#gtdb-tk-reference-data), and says which versions will work with what. Then let's say I want release95, then I go on ananconda [here](https://anaconda.org/bioconda/gtdbtk/files), and use the drop down for the Version. I try to use the highest number that will work. For example, I'm saying this because I once tried to install v.1.4.2 but only 1.4.1 was available. Anyway, once you figure out which versions exist on anaconda, create you environment [like so](https://ecogenomics.github.io/GTDBTk/installing/bioconda.html):, and make sure to name the env something informative:
    
    ```
    conda create -n gtdbtk-1.7.0 -c conda-forge -c bioconda gtdbtk=1.7.0
    
    ```
    
    Then you need to link the correct database to your program (e.g don't link r95 to v1.7.0, it needs release202):
    
    ```
    conda activate gtdbtk-1.7.0
    which gtdbtk
    echo "export GTDBTK_DATA_PATH=/path/to/release/package/" > /miniconda3/envs/gtdbtk-1.3.0/etc/conda/activate.d/gtdbtk.sh # But replace for the your path
    source /miniconda3/envs/gtdbtk-1.3.0/etc/conda/activate.d/gtdbtk.sh
    gtdbtk -h
    
    ```
    
- usearch
    
     Download the [linux binary](https://drive5.com/usearch/download.html), run `chmod +x [filename]` on that file, rename it `mv usearchvXXXX usearch`, then move the file to `/slowdata/bin/`. usearch will now be accessible anywhere. We need usearch for DasTool to work.
    
    DASTool: [https://github.com/cmks/DAS_Tool](https://github.com/cmks/DAS_Tool)
    
- DRep: https://github.com/MrOlm/drep
    
    
    Create a conda environment for dRep ([https://drep.readthedocs.io/en/latest/installation.html](https://drep.readthedocs.io/en/latest/installation.html)), then run `dRep check_dependencies` to see if the dependencies are installed. In my case, `CheckM`, `ANIcalculator` and `centrifuge` were not necessarily installed correctly. In the conda env (it should say (dRep) in front of your prompt line), type `conda install -c bioconda checkm-genome`, then rerun `dRep check_dependencies`. Repeat for ANIcalculator and centrifuge until `dRep check_dependencies` shows that all dependencies are correctly installed. PS: Don't forget to run `checkm setRoot [path to the checkm database folder]` before running DRep. ANI calculator can be downloaded here: `wget <https://ani.jgi.doe.gov/download_files/ANIcalculator_v1.tgz`>. For centrifuge follow these instructions ([https://ccb.jhu.edu/software/centrifuge/manual.shtml#obtaining-centrifuge](https://ccb.jhu.edu/software/centrifuge/manual.shtml#obtaining-centrifuge)). `wget <https://github.com/DaehwanKimLab/centrifuge/archive/refs/tags/v1.0.4.tar.gz`>
    

<aside>
😎 **We’re ready to process some data!** The methods are detailed in my paper, but here, I included the commands (”tutorial” style) in case someone would like to do these analyses on their own data.

</aside>

<aside>
⚠️ Ok it might be easier if there was a easy way to create a pipeline that would do all these processing steps! It would definitely help, but I find that each data set is unique and having control over the settings at different “steps” can be helpful. 
For example, what’s your research question? Do you care about overall taxonomic diversity or strain-level variation within a specific genus? How many different ways of binning do you need to try until it doesn’t affect your results anymore? What does your assembly quality look like (and how do you know?)? Do you want to try a couple of kmers settings to see how that’s going to affect your bins?

</aside>

# Processing the Total Metagenomes: Getting Bacterial & Archaea MAGs

![SuppFig5-Workflow-JGI-ours-combining.jpg](Lake%20Mendota%20Paper%204ca5288191e64edfab21674f31775175/SuppFig5-Workflow-JGI-ours-combining.jpg)

## Step 0. Getting a good sequence assembly

When we get raw `fastq` data back from the sequencer, we need to follow the general steps:

1. Evaluate the input reads
    1. `FastqQC` takes in reads in `R1.fastq` and `R2.fastq` format
    2. `BBTools` read correction using bbcms V. 38.86. For this dataset it was done using these commands:
    
    ```bash
    bbcms.sh -Xmx100g metadatafile=counts.metadata.json
    mincount=2 highcountfraction=0.6 in=bbcms.input.fastq.gz out1=input.corr.left.fastq.gz
    out2=input.corr.right.fastq.gz
    ```
    
2. Assemble the reads
    1. `metaSPADES` (l[ink](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5411777/#:~:text=metaSPAdes%20addresses%20various%20challenges%20of,and%20highly%20polymorphic%20diploid%20genomes.)) For this dataset, it was using these settings:
    
    ```bash
    spades.py -m
    2000 --tmp-dir cromwell_root -o spades3 --only-assembler -k 33,55,77,99,127 --meta -t 16
    -1 input.corr.left.fastq.gz -2 input.corr.right.fastq.gz.
    ```
    
3. Determine how many reads map back to the assembled contigs
    1. Can be done using any mapping tools such as `BBmap` . For this dataset, it was done using :
    
    ```bash
    bbmap.sh build=1 overwrite=true fastareadlen=500 threads=16 nodisk=true
    interleaved=true ambiguous=random rgid=filename in=reads.fastq.gz ref=reference.fasta
    out=pairedMapped.bam
    ```
    
4. Gather statistics on the assembles contigs

There are a few tutorials online, for example [here](https://training.galaxyproject.org/training-material/topics/assembly/tutorials/general-introduction/tutorial.html).

<aside>
💡 **What makes a good assembly?**
It also depends on the dataset and what we’re comparing the data to but in general you can look for:
****- **L50**: The count of smallest number of contigs whose length sum makes up half of genome size.count of smallest number of contigs whose length sum makes up half of genome size.
- **N50**: the *N50*
 is defined as the sequence length of the shortest contig at 50% of the total assembly length
- **Percentage of reads mapped back into the contigs:** A higher percentage is better (e.g. >80%)
- The **depth of coverage** of the reads mapped back onto the contigs, where ******************even coverage across the scaffold is ideal******************

</aside>

## Step 1. Binning

We will use `Metawrap` v.1.3.2 (Last updated April 20, 2021\*) to bin the assembled scaffolds into bins, and use `metabat1`, `metabat2`, and `maxbin2`.

To install metawrap, which uses python2.7, I followed the very exact instruction on the [Github page](https://github.com/bxlab/metaWRAP#installation) for the metawrap-env environment.

<aside>
💡  You can use `nohup` within a conda environment!

</aside>

The command is:

```bash
ssh patricia@sulfur.doit.wisc.edu
[enter password]

cd /slowdata/data2/Lake_Mendota_Patricia/Bacteria_Workflow/Binning
conda activate metawrap-env

for file in /slowdata/Reads/Paired_Bact_Vir_Mendota/Metagenomes/assembled_contigs/*.fna; do echo "Starting to bin $file"; metawrap binning -a $file -o ${file##*/}_with_diff_coverage -t 25 --metabat2 --metabat1 --maxbin2 --interleaved /slowdata/Reads/Paired_Bact_Vir_Mendota/Metagenomes/*.fastq; echo "Done binning $file with differential coverage, going to the next metagenome"; done

```

You can run this command for all of your 16 bacterial samples, the script is:

```bash
nohup /storage1/data10/Step1_RunMetawrap.sh > Step1_RunMetawrap.log &

```

To give you an idea of time estimate, using `25 threads`, on a `1GB assembly`, and `50GB reads file (x 30 samples for differential coverage)`, the overall run time was 24 hours. Both metabat1 abd metabat2 takes less time than maxbin2 to run. So for 16 samples, plan for at least 2 weeks without interruptions.

This version uses updated NCBI_nt, NCBI_tax and KRAKEN2 databases which were downloaded onto the server on April 20, 2021 as well.

<aside>
💡 **How does binning affect the MAGs we get?**

There are several methods for binning with all have their pros and cons, and some methods may be better suited for given datasets of various community composition complexity. 
For my purposes, I've tested binning with 
1) each metagenome using its assembly only 
2) each metagenome using the assemblies from the same environment (e.g. oxic, anoxic samples), 
3) each metagenome using all the assemblies. 

Binning with differential coverage from all the samples takes longer but results in ************more MAGs************ and of **********higher genome quality.**********

[SuppFig6--BinningComparison.pdf](Lake%20Mendota%20Paper%204ca5288191e64edfab21674f31775175/SuppFig6--BinningComparison.pdf)

</aside>

## Step 2. Bin refinement

Instead of using the `bin_refinement` available through the `metawrap` program, we will use [DASTool](https://github.com/cmks/DAS_Tool). For folder organization, I like having a DASTool Folder within each of my `metawrap` folder,  but another method could be to keep the analyses completely separate.

Before running `DAS_Tool`, we need to make table that tells the program the links between the scaffolds ID and the bins they come from. Make sure that `usearch` is functional.

<aside>
💡 Note that in this particular scenario, we also had a separate binning done by JGI (metabat2 with reads from that single assembly) I want to refine using all the bin sets, that means MAGs processed by them too, and not just those by `metawrap`.

</aside>

An example command is:

```bash
cd /storage1/data10/Metagenome/Ga0485158_contigs.fna_with_diff_coverage/
#Create a folder called DAS_Tool INSIDE of the binning folders.
mkdir DAS_Tool
cd DAS_Tool
#Add the JGI mags
mkdir ../jgi_mags
#Copy over the JGI MAGS
cp /storage1/data10/Metagenome/JGI_MAGS_normal/Ga0485158/* ../jgi_mags/.
cp /storage1/data10/Metagenome/JGI_MAGS_CPR/Ga0485158/* ../jgi_mags/.
#Use Fastq to Scaffolds to make a TSV file with 2 columns pairing the scaffolds and the bin names.
/slowdata/archive/DAS_Tool-master/src/Fasta_to_Scaffolds2Bin.sh -e fa -i ../maxbin2_bins/ > maxbin2_scaf2bin.tsv
/slowdata/archive/DAS_Tool-master/src/Fasta_to_Scaffolds2Bin.sh -e fa -i ../metabat1_bins/ > metabat1_scaf2bin.tsv
/slowdata/archive/DAS_Tool-master/src/Fasta_to_Scaffolds2Bin.sh -e fa -i ../metabat2_bins/ > metabat2_scaf2bin.tsv
/slowdata/archive/DAS_Tool-master/src/Fasta_to_Scaffolds2Bin.sh -e fna -i ../jgi_mags/ > jgi_bins_scaf2bin.tsv
#Run DAS_Tool:
DAS_Tool --write_bins 1 --write_bin_evals 1 --create_plots 1 -i metabat1_scaf2bin.tsv,metabat2_scaf2bin.tsv,maxbin2_scaf2bin.tsv,jgi_bins_scaf2bin.tsv -l metabat1,metabat2_ours,maxbin,metabat2_jgi -c /slowdata/Reads/Paired_Bact_Vir_Mendota/Metagenomes/assembled_contigs/Ga0485158_contigs.fna -o DASTool_output_with_jgi --threads 15
cd  ..
```

Repeat this command for you 15 other remaining folders (16 in total).

Note: I have a script that will do that:
`nohup bash /storage1/data10/Metagenome/Step2_RunDastool.sh > Step2_RunDastool.log &`

After DAS_Tool has finished running type this to view how many bins were created for each sample after refinement:

```bash
grep "above score threshold" /storage1/data10/Metagenome/**/DAS_Tool/*.log
```

## Step 3. Spitting up CPR from non CPR bins

Completeness and contamination of the MAGs can be calculated using [CheckM](https://github.com/Ecogenomics/CheckM). Here, instead of running the CheckM command 16 times, I make a new folder named "DASTool_refined_bins_all" and I can check all my MAGS at the same time there.

Example:

```bash
cd /storage1/data10/Metagenome/Ga0485159_contigs.fna_with_diff_coverage/DAS_Tool/DASTool_output_DASTool_bins/
for f in *.fa; do mv "$f" "Ga0485159_${f%.fa}.fa"; done
```

**This is very important** because I don't want to mix up where the MAGs are coming from. Also don't forget to set the format of your table to **tab separated** (`--tab_table`)! It will be easier to parse later.

I have a script for that:
`Step3.1_Rename_refined_bins.sh` and
`Step3.1.1_copy_refined_renamed_bins_to_new_folder.sh`

Once the bins are renamed and moved to a folder (e.g. `Refine_bins_including_jgi_bins/.`), we must run GTDB-tk first. 

**Why?** Because we will put the CPR bins together because dereplicating, and the non-CPR bins together before dereplicating. Since `dRep` will rely on `checkM`, and since `checkM`  has 2 sets of markers genes, we **cannot** just use 1 checkM run on ***all*** the MAGS, since the completeness values won't be as accurate between the CPR and non-CPR.

## Step 4. Taxonomic assignment using GTDB-tk

To run GTDBTK on *all the genomes*:

```bash
gtdbtk classify_wf --genome_dir Refined_bins_including_jgi_bins -x fa --cpus 15 --out_dir gtdbtk_output_on_refined_bins_after_incl_jgi_bins
```

To get the genomes with p__Patescibacteria in the name:

```bash

grep 'p__Patescibacteria' storage1/data10/Metagenome/gtdbtk_output_on_refined_bins_after_incl_jgi_bins//storage1/data10/Metagenome/gtdbtk_output_on_refined_bins_after_incl_jgi_bins/gtdbtk.bac120.summary.tsv
```

Then we move the CPR bins in a folder and the non CPR bins in another folder:
`Step4_run_GTDBTK_on_Refined_bins_including_jgi_bins.sh`

## Step 5. Splitting up CPR and non-CPR bins:

Now that we know which bins are CPR vs non CPR let's copy them to new folders:
`Step5.1_copy_cpr_bins_after_new_refinement.sh` and `Step5.1.1_copy_normal_bins_after_new_refinement.sh`

## Step 6. Quality checking on CPR vs non CPR bins:

After activating the checkm environment (i.e. `conda activate CheckM`), run CheckM normally (`lineage_wf`) on the Refined bin set, then use the `analyze` function to use the `cpr_check.hmm` marker genes on it, and then use `qa` to get the tab delimited table:

```bash
checkm lineage_wf ./Refined_bins_including_jgi_bins_CPR/ ./Checkm_on_Refined_Bins_CPR_new -t 15 --pplacer_threads 15 --tmpdir /storage1/data10/tmp --tab_table -f checkm_on_refined_bins_CPR_new.tsv -x fa

checkm analyze /slowdata/databases/cpr_checkm.hmm -x fa -t 15 --tmpdir /storage1/data10/tmp/ /storage1/data10/Metagenome/Refined_bins_including_jgi_bins_CPR/ /storage1/data10/Metagenome/Checkm_on_Refined_Bins_CPR_new/

checkm qa /slowdata/databases/cpr_checkm.hmm --tab_table -f checkm_on_refined_bins_CPR_marker_set_new.tsv -t 15 --tmpdir /storage1/data10/tmp/ /storage1/data10/Metagenome/Checkm_on_Refined_Bins_CPR_new/

```

The table of interest is `checkm_on_refined_bins_CPR_new.tsv`

The script is `Step6.1_checkm_on_refined_bins_CPR_new.sh`

## Step 7. Dereplicate the CPR bin set:

<aside>
💡 **There are some essential steps to enable using the checkM TSV file in dRep that is not totally obvious:**
You will have to do these modifications:

</aside>

<aside>
💡 Replace headers by `genome`, `completeness`, and `contamination` instead of them being capitalized.

</aside>

<aside>
💡 Make it a `CSV` file instead of a `TSV` file

</aside>

<aside>
💡 Add the file extension (.fasta, .fna, etc) to the genome names in the first column/

</aside>

After running these modification, then you can run dRep:
The argument `---genomeInfo` is where you provide it the TSV (actually a CSV) file.

```bash
dRep dereplicate -p 10 -comp 50 --genomeInfo checkm_on_refined_bins_CPR_marker_set_new.tsv output_dir_dRep_new -g ./Refined_bins_including_jgi_bins_CPR/*.fa

```

Script: `Step7.1_dRep_on_CPR_bins_new.sh`

## Step 8 & 9. Repeat for the non-CPR bins:

Similarly, we will run `CheckM` and `dRep` on the non-CPR bins sets. No need to run GTDBTK again just yet, because we ran GTDBTK on the full set of refined genomes already. Just make sure that after moving the CPR bins outside of that folder, you have a folder with non-CPR bins only.

```bash
checkm lineage_wf ./Refined_bins_including_jgi_bins_normal/ ./Checkm_on_refined_bins_normal_new -t 25 --pplacer_threads 25 --tmpdir /storage1/data10/tmp --tab_table -f checkm_on_refined_bins_split_new.tsv  -x fa

```

Now **do the modifications to the CheckM TSV file**, and then:

```bash
dRep dereplicate -p 20 -comp 50 -con 10 --genomeInfo checkm_on_refined_bins_split_new.tsv output_dir_dRep_normal_new -g ./Refined_bins_including_jgi_bins_normal/*.fa
```

Scripts are here: `Step8.1_checkm_on_refined_bins_new.sh` and `Step9.1_dRep_on_refined_bins_new.sh`

## Step 10. Final Taxonomic Assignment

After refinement and dereplication, some bin names might have changed (notice names ending with "sub")

Let's run GTDTK one more time

Copy over all your final bins in `Final_Bin_Set/`

then:

```bash
gtdbtk classify_wf --genome_dir Final_Bin_Set/ --out_dir GTDBTK_on_Final_Bin_Set -x fa --cpus 25 --pplacer_cpus 25  --tmpdir /storage1/data10/tmp/

```

Which is Script : `Step10_GTDBTK_on_Final_Bin_Set.sh`

# 🥳 The basic processing of MAGs is done! 🥳

At this point, you have a collection of MAGs and depending on the analyses and questions, you can do different things.