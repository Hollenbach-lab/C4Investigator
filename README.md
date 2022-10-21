# C4Investigator
An R-based bioinformatic pipeline to determine complement component 4 (C4) copy number and genotypes from short-read genomic sequencing data.


## Language
R


## System compatibility
* Linux (tested on Ubuntu, CentOS)
* OS X
* Windows (untested)


## System dependencies download and install
* bowtie2 (tested with v2.4.1, but should be compatible with newer versions)
* samtools (tested with v1.11, but should be compatible with newer versions)


### 1. Download C4Investigator and run a 1000 Genomes population

```shell
git clone https://github.com/wesleymarin/C4Investigator.git --single-branch --branch master
```

### 2. Download 1000Genomes data to test

```shell
cd C4Investigator/ && \
Rscript 1000Genomes_download_coordinator.R --population FIN --fqDirectory <fqDirectory> --threads <threads>
```
`--fqDirectory`     Directory to download FQ files to

`--threads`         Number of compute threads to utilize


### 3. Run on 1000Genomes data, or your own data (starting with fastq files)

```shell
cd C4Investigator/ &&
Rscript C4Investigator_run.R --fqDirectory <FQ_directory> --resultsDirectory <OutputDirectory> --fastqPattern <fq/fastq> --threads <threads>
```

`--fqDirectory`     Directory holding your fastq data.

`--resultsDirectory`  Desired output directory

`--fastqPattern`      A string that is shared across all of your fastq file names (used to find fq files and match pairs), this is usually fq or fastq

`--threads`          Number of compute threads to utilize


### 4. Output
`C4Investigator_c4_summary.csv`   <- this is the main output file that contains C4 copy number calls.

`C4Investigator_c4_detailed.csv`  <- This file contains more detailed information that is used to determine C4 copy number calls.

`pileups/C4_DP_<sampleID>.csv`    <- This is the alignment pileup file used for SNP calling
                                
`pileups/C4_SNP_<sampleID>.csv`   <- This is the file containing C4 SNP calls
                                 
`plots/<sampleID>_c4_plot.html`   <- A web viewable, interactive plot of the C4 alignment
                                 
                                 
