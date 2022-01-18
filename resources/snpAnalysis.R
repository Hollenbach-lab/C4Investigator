library(data.table)
library(UniprotR)
library(stringr)
library(seqinr)
library(Rsamtools)
library(GenomicAlignments)
library(gtools)

setwd('E:/workData/C4Investigator/')
source('c4_analysisSupport_functions.R')

output.dir <- 'output/tgpData/'
rawData.dir <- 'tgpData'

# Resources ---------------------------------------------------------------

#C4B.uniprot <- processUniprot('P0C0L5')
#C4A.uniprot <- processUniprot('P0C0L4')

mrna.manifest.dt <- fread('tgpData/mRNA_manifest.tsv',sep='\t')
manifest.dt <- fread('tgpData/manifest.tsv',sep='\t')
fullData.dt <- fread('output/tgpData/fullData.csv',header=T)
fullData.dt <- fullData.dt[,-1]

hom.snpDT <- fread('output/tgpData/homSampleSNP.csv',header = T)
hom.snpDT <- hom.snpDT[,-1]
het.snpDT <- fread('output/tgpData/hetSampleSNP.csv',header = T)
het.snpDT <- het.snpDT[,-1]

exonPos.vect <- grep( 'E', setdiff( colnames(hom.snpDT), 'sampleID' ), value=T)
protLen <- length(exonPos.vect)/3

codonConv.list <- list()
initPos <- 1
for( codonPos in paste0('C_',1:protLen) ){
  posVect <- initPos:(initPos+2)
  codonConv.list[[codonPos]] <- exonPos.vect[posVect]
  initPos <- initPos+3
}

hetPos.vect <- names( which( apply(hom.snpDT[,exonPos.vect,with=F],2,num_unique_nuc) > 1 ) )



# C4A sequence analysis ---------------------------------------
homIDvect.c4a <- fullData.dt[C4B_copy == 0 & C4A_copy > 0]$.id

temp.homVect.c4a <- c()
temp.hetVect.c4a <- list()
for( homID in homIDvect.c4a ){
  sample.snpDT <- hom.snpDT[sampleID == homID,1:(ncol(hom.snpDT)-1)]
  
  hetPos.boolVect <- apply(sample.snpDT, 2, num_unique_nuc) > 1
  hetPos.int <- sum(hetPos.boolVect)
  
  if( hetPos.int == 0){
    temp.homVect.c4a <- c(temp.homVect.c4a, homID)
  }else{
    temp.hetVect.c4a[[homID]] <- list('hetPos.N'=hetPos.int,'hetPos.vect'=names(which(hetPos.boolVect)))
  }
  
  cat('\n',homID,'',hetPos.int)
}


table(unlist(sapply(temp.hetVect.c4a, function(x)x$hetPos.vect)))


for( homID in temp.homVect.c4a ){
  cat('\n\n',homID)
  sample.snpDT <- hom.snpDT[sampleID == homID,1:(ncol(hom.snpDT)-1)]
  sample.snpDT <- sample.snpDT[1,]
  
  nuc.vect <- as.character(sample.snpDT)
  protStr <- gsub('*','',paste0(translate(nuc.vect),collapse=''),fixed=T)
  
  mismatchPos.vect <- which( strsplit(protStr,'')[[1]] != strsplit(C4A.uniprot$Sequence,'')[[1]] )
  
  if( length(mismatchPos.vect) > 0 ){
    for( mismatch.pos in mismatchPos.vect ){
      
      mismatchCodon.list <- codonConv.list[mismatch.pos]
      
      sample.mismatchPos.nuc <- strsplit(protStr,'')[[1]][mismatch.pos]
      reference.mismatchPos.nuc <- strsplit(C4A.uniprot$Sequence,'')[[1]][mismatch.pos]
      
      cat('\n\t',paste0(mismatch.pos,': ',reference.mismatchPos.nuc,' -> ',sample.mismatchPos.nuc))
    }
  }
}




# C4B sequence analysis ---------------------------------------------------
homIDvect.c4b <- fullData.dt[C4A_copy == 0 & C4B_copy > 0]$.id

temp.homVect.c4b <- c()
temp.hetVect.c4b <- c()
for( homID in homIDvect.c4b ){
  sample.snpDT <- hom.snpDT[sampleID == homID,1:(ncol(hom.snpDT)-1)]
  
  hetPos.int <- sum(apply(sample.snpDT, 2, num_unique_nuc) > 1)
  
  if( hetPos.int == 0){
    temp.homVect.c4b <- c(temp.homVect.c4b, homID)
  }else{
    temp.hetVect.c4b <- c(temp.hetVect.c4b, homID)
  }
  
  cat('\n',homID,'',hetPos.int)
}

for( homID in temp.homVect.c4b ){
  cat('\n\n',homID)
  sample.snpDT <- hom.snpDT[sampleID == homID,1:(ncol(hom.snpDT)-1)]
  sample.snpDT <- sample.snpDT[1,]
  
  nuc.vect <- as.character(sample.snpDT)
  protStr <- gsub('*','',paste0(translate(nuc.vect),collapse=''),fixed=T)
  
  mismatchPos.vect <- which( strsplit(protStr,'')[[1]] != strsplit(C4B.uniprot$Sequence,'')[[1]] )
  
  if( length(mismatchPos.vect) > 0 ){
    for( mismatch.pos in mismatchPos.vect ){
      
      mismatchCodon.list <- codonConv.list[mismatch.pos]
      
      sample.mismatchPos.nuc <- strsplit(protStr,'')[[1]][mismatch.pos]
      reference.mismatchPos.nuc <- strsplit(C4B.uniprot$Sequence,'')[[1]][mismatch.pos]
      
      cat('\n\t',paste0(mismatch.pos,': ',reference.mismatchPos.nuc,' -> ',sample.mismatchPos.nuc))
    }
  }
}


'
Rare C4B del variant, described https://www.frontiersin.org/articles/10.3389/fimmu.2021.739430/full
E17_192.C deletion
  not detected in this dataset.

C4B del variant
E15_110 G -> A (premature stop codon)
  2 samples positive
  "HG00626" "NA18618"
  Both EAS population, consistent with Frontiers paper
'




View(hom.snpDT[sampleID %in% homIDvect.c4a])


# Phasing testing ---------------------------------------------------------

'
Want to add code to C4 investigator to test alternative normalizations
  Getting average depth of Chr6 / WGS
    samtools coverage -r chr6 http://s3.amazonaws.com/1000genomes/1000G_2504_high_coverage/data/ERR3239597/NA19066.final.cram -H -o cov.dump
  
  Build bash command for getting average depth from bam file
  
Want to save BAM files for paired-end read phasing
'

resourcesDirectory <- 'resources/'
## Setting up C4 reference resources
referencePath <- file.path(resourcesDirectory,'all_onelines_oneDel_bShort') ## C4 alignment reference
referenceKeyPath <- file.path(resourcesDirectory,'reference_key.txt') ## Key for interpreting sequence names in the reference file
alignedC4Path <- file.path(resourcesDirectory,'c4only_onelines_oneDel_bShort.fasta') ## C4Along and C4Bshort alignment reference alleles matched by position
alignedC4Path <- normalizePath(alignedC4Path,mustWork=T) ## Making sure the matched reference allele file exists

## Reading in the C4 reference key
referenceKeyDF <- read.table(referenceKeyPath, stringsAsFactors = F) ## Read in the file as a table
alignedLocusVect <- c(unique(referenceKeyDF[,6]),'C4') ## Pull out the locus names, adding in overall C4

## Build up a conversion list to convert genome locus names to common gene names
referenceKeyList <- list() 
for(i in 1:nrow(referenceKeyDF)){
  referenceKeyList[[referenceKeyDF[i,1]]] <- referenceKeyDF[i,6]
}

## Read in the aligned C4 fasta (a dataframe for each locus)
c4AlleleDF <- read.c4_dataframe_from_reference_fasta(alignedC4Path, referenceKeyList)

c4.charList <- c4AlleleDF[1,,drop=T]
c4.exonCoordVect <- names(c4.charList)[unlist(lapply(c4.charList, general.identify_uppercase))]

## Build up deletion index lists, which are needed for converting read coordinates between Along and Bshort
inverseDeletionIndexList <- build.c4_inverse_deletion_index_list(c4AlleleDF)
deletionIndexList <- build.c4_deletion_index_list(c4AlleleDF)






# Processing BAM file to generate phase matrix ----------------------------
currentSample.id <- 'HG01345'
toPhase.samDT <- c4.process_samFile_to_phaseDT( currentSample.id, rawData.dir )

currentSample.hetPosVect <- c(AB.snpPosVect,RgCh.snpPosVect)
targetRegion.posVect <- names( c4.feature.vect[c4.feature.vect['E25_1']:c4.feature.vect['I28_20']] )
targetRegion.hetPosVect <- intersect(currentSample.hetPosVect, targetRegion.posVect)

toPhase.snpVect <- c4.feature.vect[targetRegion.hetPosVect]
if( currentSample.id %in% hom.snpDT$sampleID ){
  toPhase.hetSnpDT <- hom.snpDT[ sampleID == currentSample.id ][,names(toPhase.snpVect),with=F]
}else{
  toPhase.hetSnpDT <- het.snpDT[ sampleID == currentSample.id ][,names(toPhase.snpVect),with=F]
}

phase.mat <- c4.samDT_to_phaseMAT(toPhase.samDT, toPhase.snpVect, toPhase.hetSnpDT )

AB.RgCh.inPhase.bool <- c4.verify_AB_RgCh_phase( phase.mat )

if(!AB.RgCh.inPhase.bool){stop('AB / RgCh not in phase')}

# Processing phased SNPs --------------------------------------------------
phased.tab <- c4.phase_AB_RgCh( phase.mat )

phased.tab






















hetPos <- 'E3_157'
# Hom sample analysis
for(hetPos in hetPos.vect){
  snpCount.tab <- table( hom.snpDT[,hetPos,with=F] )
  snpCount.tab <- snpCount.tab[ is_nuc(names(snpCount.tab)) ]
  
  if( sum( snpCount.tab > 1 ) > 1 ){
    cat('\n',hetPos,'\t',names(snpCount.tab),'\t',snpCount.tab)
    codonPos <- names( which( sapply(codonConv.list, function(x) hetPos %in% x) ) )
    codonVect <- codonConv.list[[codonPos]]
  }
}
hom.snpDT[,codonVect,with=F]

# Het sample analysis
for(hetPos in hetPos.vect){
  snpCount.tab <- table( het.snpDT[,hetPos,with=F] )
  snpCount.tab <- snpCount.tab[ is_nuc(names(snpCount.tab)) ]
  
  if( sum( snpCount.tab > 1 ) > 1 ){
    cat('\n',hetPos,'\t',names(snpCount.tab),'\t',snpCount.tab)
    codonPos <- names( which( sapply(codonConv.list, function(x) hetPos %in% x) ) )
    codonVect <- codonConv.list[[codonPos]]
  }
}

'
E36_4
  CCC>TCC Variant not in NCBI
  P>S | NP > P
  13/95 in hom dataset
    13 are AFR
  82/2902 in het dataset
    78 are AFR
    4 are AMR
  samples with mRNA
     "NA18489" "NA18908" "NA18508" "NA19223" "NA19146" "NA18923" "NA19096" "NA19116" "NA18909" "NA19185" "NA19257" "NA18486" "NA19092" "NA19206" "NA19131"
     
E33_38
  GAT > GAC
  D = D

E33_6
  GCA > CCA
  A > P | NP > NP
  15/95 in hom dataset
    15 are AFR
  94/2902 in het dataset
    91 are AFR
    3 are AMR
  samples with mRNA
    "NA18489" "NA18508" "NA19189" "NA19223" "NA19146" "NA19096" "NA19116" "NA18909" "NA19185" "NA19257" "NA18486" "NA19206" "NA19131"

E29_180
  GCG > TCG | NCBI: rs9501603
  A > S | NP > P
  11/95 in hom dataset
    2 EAS
    2 EUR
    7 SAS
  331/2902 in het dataset
    28 AFR
    37 AMR
    78 EAS
    118 EUR
    80 SAS
  mRNA
    HG00118 HG00120 HG00149 HG00151 HG00238 HG00097 HG00105 HG00310 HG00143 HG00181 HG00232 HG00263 HG00268 HG00275 HG00282 NA11832 NA11995 NA12044 NA11831 NA11881 NA11918 NA11994 NA12005 NA12889 NA12873 NA20585 NA20807 NA20515 NA20527 NA20539 NA20541 NA20798 NA20806 NA20813 HG00157 HG00253 HG00277 HG00107 HG00142 HG00114 HG00121 HG00159 HG00178 HG00145 HG00236 HG00332 HG00337 HG00344 HG00361 HG00351 HG00382 HG00250 HG00313 HG01791 HG00189 HG00379 NA12413 NA12829 NA12843 NA06985 NA12275 NA12340 NA12383 NA11892 NA11931 NA12004 NA07048 NA19248 NA12763 NA12814 NA19149 NA20759 NA20778 NA20785 NA12717 NA12750 NA11930 NA20518 NA20796 NA20809

E29_59
  CCG > CCA
  P = P

E29_50
  TCG > TCA
  S = S

E28_23
  AGC > AAC
  S > N | P > P
  Very common

E26_55
  GGC > GGA
  G = G

E25_64
  GAC > GGC | rs2258218 & rs147162052
    rs2258218: no publications
    rs147162052: no publications
    
  D > G | A > NP
  Very common

E24_12
  TTG > CTG
  L = L

E23_23
  GAA > GAC
  E > D | A > A
  
  
E21_127
  GCC > ACC | rs796750528 & rs429329
    rs796750528: no publications
    rs429329: no publications
    
  A > T | NP > P
  
  19/95 hom samples
    AFR EAS EUR SAS 
     4   2   6   7 
     
  424/2902 het samples
    AFR AMR EAS EUR SAS 
    57  54  59 182  72 
    
  mRNA samples
    HG00327 HG00106 HG00118 HG00120 HG00137 HG00149 HG00151 HG00238 HG00257 HG00097 HG00100 HG00105 HG00117 HG00310 HG00143 HG00101 HG00179 HG00232 HG00345 HG00251 HG00263 HG00268 HG00275 HG00282 HG02215 NA11832 NA11894 NA11995 NA12044 NA11831 NA12342 NA11881 NA11918 NA11994 NA12889 NA12842 NA18488 NA12873 NA12872 NA07051 NA07056 NA18868 NA19114 NA20585 NA20802 NA20807 NA20515 NA20527 NA20539 NA20541 NA20798 NA20806 NA20813 HG00157 HG00253 HG00265 HG00272 HG00277 HG00284 HG00109 HG00099 HG00107 HG00142 HG00114 HG00121 HG00159 HG00178 HG00145 HG00243 HG00332 HG00337 HG00344 HG00349 HG00351 HG00356 HG00380 HG00382 HG00250 HG00306 HG00313 HG00096 HG01791 HG00158 HG00189 HG00266 HG00280 HG00146 HG00336 NA12413 NA12761 NA12829 NA12843 NA06985 NA12155 NA12275 NA12282 NA12340 NA11892 NA11931 NA12004 NA18923 NA07000 NA19236 NA20536 NA20543 NA20586 NA12751 NA12763 NA12814 NA20783 NA19149 NA18867 NA20588 NA20759 NA07357 NA20773 NA20778 NA20785 NA20800 NA12154 NA12717 NA12750 NA11930 NA18873 NA19117 NA20518 NA20758 NA20796 NA20809 NA20506
  
  
E20_23
  GTT > GTC
  V = V
  
E17_106
  CCG > CTG
  P > L | NP = NP

E13_122 | C_549
  CAT > CCT
    rs2229405: http://hdl.handle.net/20.500.11937/1013
    
  H > P | B > NP
  4/94 hom samples
    AFR AMR EAS 
     1   1   2 
  518/2902 het samples
    AFR AMR EAS EUR SAS 
     14  81 244  89  90
  HG00339 HG00341 HG00346 HG00353 HG00358 HG00118 HG00237 HG00321 HG00364 HG00371 HG00268 NA06984 NA11894 NA12006 NA11920 NA12830 NA12842 NA12827 NA19204 NA20504 NA20516 NA20542 NA20798 NA20801 HG00176 HG00234 HG00260 HG00119 HG00126 HG00309 HG00323 HG00320 HG00328 HG00335 HG00349 HG00366 HG00356 HG01789 HG00160 HG00141 HG00379 NA12348 NA12778 NA12874 NA10847 NA12249 NA07346 NA20529 NA20531 NA20757 NA12751 NA12763 NA12814 NA20540 NA20773 NA20778 NA12750 NA11992 NA20513 NA20532 NA20544 NA20506

E12_144
  GCC > GCT
  A = A

E11_117
  GAC > GAT
  D = D

E9_128
  TCT > TAT
  S > Y | P = P

E9_105
  TAC > TAT
  Y = Y

E3_157
  CTC > GTC
  L > V | NP = NP
'

names( which( sapply(codonConv.list, function(x) 'E13_122' %in% x) ) )

#hetPos <- 'E29_150'
hetNuc <- 'T'
table( manifest.dt[ manifest.dt[,'Sample name'][[1]] %in% hom.snpDT$sampleID[hom.snpDT[,hetPos,with=F] == hetNuc], ]$`Superpopulation code`)
table( manifest.dt[ manifest.dt[,'Sample name'][[1]] %in% het.snpDT$sampleID[het.snpDT[,hetPos,with=F] == hetNuc], ]$`Superpopulation code`)
cat(mrna.manifest.dt[,'Sample name'][[1]][mrna.manifest.dt[,'Sample name'][[1]] %in% c( het.snpDT$sampleID[het.snpDT[,hetPos,with=F] == hetNuc], hom.snpDT$sampleID[hom.snpDT[,hetPos,with=F] == hetNuc])])
