# C4Investigator
An R-based bioinformatic pipeline to determine complement component 4 (C4) copy number and genotypes from short-read genomic sequencing data.


## Language
R


## System compatibility
* Linux (tested on Ubuntu, CentOS)
* OS X
* Windows (untested)



### 1. Download C4Investigator 

```shell
git clone https://github.com/Hollenbach-lab/C4Investigator.git
cd C4Investigator
```

### 2. Build the Singularity container
```shell
singularity build --fakeroot c4investigator.sif c4investigator.def
```

### 3. Download 1000Genomes data to test

```shell
singularity exec c4investigator.sif Rscript 1000Genomes_download_coordinator.R --population FIN --fqDirectory <fqDirectory> --threads <threads>
```
`--fqDirectory`     Directory to download FQ files to

`--threads`         Number of compute threads to utilize


### 4. Run on 1000Genomes data, or your own data (starting with fastq files)

```shell
singularity exec c4investigator.sif Rscript C4Investigator_run.R --fqDirectory <FQ_directory> --resultsDirectory <OutputDirectory> --fastqPattern <fq/fastq> --threads <threads>
```

`--fqDirectory`     Directory holding your fastq data.

`--resultsDirectory`  Desired output directory

`--fastqPattern`      A string that is shared across all of your fastq file names (used to find fq files and match pairs), this is usually fq or fastq

`--threads`          Number of compute threads to utilize


### 5. Output
`C4Investigator_c4_summary.csv`   <- this is the main output file that contains C4 copy number calls.

`C4Investigator_c4_detailed.csv`  <- This file contains more detailed information that is used to determine C4 copy number calls.

`pileups/C4_DP_<sampleID>.csv`    <- This is the alignment pileup file used for SNP calling
                                
`pileups/C4_SNP_<sampleID>.csv`   <- This is the file containing C4 SNP calls
                                 
`plots/<sampleID>_c4_plot.html`   <- A web viewable, interactive plot of the C4 alignment
                                 
                                 
# Troubleshooting
Please send a copy of your error message to rayo.suseno@ucsf.edu, or raise an issue in this repository. We would be happy to receive community feedback as we aim to continuously improve our pipeline.