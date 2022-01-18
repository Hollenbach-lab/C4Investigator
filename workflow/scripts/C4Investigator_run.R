# ---- DEPENDENCIES ----
library(data.table)
library(methods)
library(stringr)
library(plotly)



##' Checks that the input is well-formed and returns it or errors out
assert_not_null <- function(x){
  stopifnot(!is.null(x),length(x)>0)
  return(x)
}

workingDirectory <- assert_not_null(snakemake@params[["workingDirectory"]])
threads <- as.integer(snakemake@threads)
resultsDirectory <- assert_not_null(snakemake@params[["resultsDirectory"]])
minDP <- as.numeric(assert_not_null(snakemake@params[["minDP"]]))
projectName <- assert_not_null(snakemake@params[["projectName"]])
maxReadThreshold <- as.integer(assert_not_null(snakemake@params[["maxReadThreshold"]]))
output_f <- assert_not_null(unlist(snakemake@output))
inputf1 <- normalizePath(assert_not_null(snakemake@input[["fq1"]]),mustWork=TRUE)
inputf2 <- normalizePath(assert_not_null(snakemake@input[["fq2"]]),mustWork=TRUE)

owd <- getwd()
setwd(workingDirectory)

source('resources/general_functions.R')
source('resources/c4Copy_functions.R')
source('resources/generalAlignment_functions.R')


## Creating the sample object class
sample <- setRefClass("sample",
                      fields=list(name='character',
                                  fastq1path='character',
                                  fastq2path='character',
                                  gzip='logical',
                                  samPath='character',
                                  bamPath='character',
                                  failed='logical'))


sessionInfo()
setDTthreads(threads)
# Preparation -------------------------------------------------------------
outDir <- pathObj(name='output_directory',path=resultsDirectory)
outDir$dirGen()

# Build up a list of sample objects


sampleList <- list(sample(name=assert_not_null(snakemake@params[["samplename"]]),
                          fastq1path=inputf1,
                          fastq2path=inputf2,
                          gzip=TRUE,
                          failed=FALSE))

names(sampleList) <- snakemake@params[["samplename"]]



c4.run_script <- function( sampleList, projectName, resultsDirectory, threads, maxReadThreshold, minDP ){
  
  
  ## Setting up the results directory
  dir.create(resultsDirectory,showWarnings = F) ## This function will output an inconsequential warning if the directory already exists
  resultsDirectory <- normalizePath(resultsDirectory, mustWork = T) ## Make sure the results directory exists
  
  ## Read in probe CSV
  #c4.probeDF <- read.csv(file=file.path('resources/C4_exon_frame_shift_mut_probes.csv'),check.names=F,
  #                       stringsAsFactors = F)
  #row.names(c4.probeDF) <- c4.probeDF$ProbeID
  

  
  c4.initialize_indel_file <- function(resultsDirectory){
    path <- file.path( resultsDirectory, 'indel.bySample.txt')
    textStr <- ''
    cat(textStr, file=path)
    return(path)
  }
  
  c4.initialize_phasing_file <- function(resultsDirectory){
    path <- file.path( resultsDirectory, 'snpPhasing.importantSNPs.bySample.txt')
    textStr <- ''
    cat(textStr, file=path)
    return(path)
  }
  
  #indel.path <- c4.initialize_indel_file(resultsDirectory)
  
  #phase.path <- c4.initialize_phasing_file(resultsDirectory)
  
  ## Initialize a list for storing C4 alignments
  c4BuildList <- list()
  
  ## Initialize a dataframe for storing various depth values
  resultsDF <- data.frame(matrix(0,nrow=length(sampleList),ncol=18),row.names = names(sampleList),stringsAsFactors=F)
  colnames(resultsDF) <- c('median_c4mid','median_c4del','mean_c4ins','median_tnxb','mean_c4a_g1', 'mean_c4a_g2','mean_c4b_g1', 'mean_c4b_g2',
                           'mean_c4a', 'mean_c4b', 'c4aL', 'c4bL', 'median_mhc_c4', 'median_mhc_tnxb', 'mhc_c4_copy','c4_exon_del','c4_exon29_ins','c4_exon29_norm')
  
  #cffCountDF <- data.frame(matrix(0,nrow=length(sampleList),ncol=nrow(c4.probeDF)), row.names=names(sampleList), stringsAsFactors = F)
  #colnames(cffCountDF) <- rownames( c4.probeDF )
  
  for(currentSample in sampleList){
    currentSample <- bowtie2.c4_alignment(bowtie2, referencePath, threads, currentSample, resultsDirectory)
      
    currentSample <- samtools.sam_to_bam(samtools, currentSample, resultsDirectory, threads)
      
    currentSample <- bowtie2.mhc.c4_alignment(bowtie2, mhcReferencePath, threads, currentSample, resultsDirectory)
      
    currentSample <- samtools.mhc.sam_to_bam(samtools, currentSample, resultsDirectory, threads)
    
    # BAM sort
    currentSample <- samtools.mhc.sort(samtools, currentSample, resultsDirectory, threads)
    # depth
    currentSample <- samtools.all_depth(samtools, currentSample, resultsDirectory)
    
    # plot
    depthTable <- tryCatch({
      read.table(currentSample$allDepthPath, header=F, stringsAsFactors = FALSE)
    },error=function(x){
      return(character(0))
    }
    )
    
    if(length(depthTable) == 0 ){
      cat('\nBad data, skipping.')
      next
    }
    
    depthTable$V2 <- 28510119 + depthTable$V2 
    
    c4 <- depthTable[((depthTable$V2 <32035277) & (depthTable$V2 > 31982108)),]
    TNXB <- depthTable[((depthTable$V2 <32098198) & (depthTable$V2 > 32041349)),]
    c4Median <- median(c4$V3)
    TNXBMedian <- median(TNXB$V3)
    
    if( c4Median < minDP ){
      cat('\nBad data, skipping')
      next
    }
    
    # /s -----------------------------------------------------------
    
    
    
    ## Count how many header lines there in the SAM file so they can be skipped during read in
    headerLineCountInt <- samfile.count_header_lines(currentSample)
    
    ## Read in the SAM file as a data table
    samTable <- samfile.read_whole_genome_sam(currentSample$samPath, headerLineCountInt, referenceKeyList)
    
    file.remove(currentSample$samPath) ## Remove the SAM file to save space
    file.remove(currentSample$mhcSamPath) ## Remove the SAM file to save space
    
    
    
    ## Count CFF matches
    #c4.mutCountList <- c4.count_cff_matches(currentSample, c4.probeDF, samTable)
    #cffCountDF[currentSample$name,names(c4.mutCountList)] <- c4.mutCountList
    #write.csv(cffCountDF, file=file.path(resultsDirectory,paste0(projectName,'_c4_cffCount.csv')))
    

    
    #samTable$reference_name <- tstrsplit(samTable$reference_name, '_')[[3]]
    
    ## Use the reference name conversion list to apply common gene names
    #samTable$locus <- unlist(lapply(samTable$reference_name, function(x){referenceKeyList[[x]]}))
    
    ## Pull out an index of reads that aligned to C4
    #C4ReadsPosVect <- grep('C4',samTable$locus,fixed=T)
    
    ## Initialize a column of the data table for storing whether the read aligned to C4 or not
    #samTable$isC4 <- F
    #samTable$isC4[C4ReadsPosVect] <- T ## Set the C4 bool to T for all C4 indexed reads
    
    ## Pull out the subset of C4 aligned reads as its own data table
    c4SamTable <- samTable[isC4 == T]
    
    ## Intialize and set a column for read lengths
    c4SamTable$readLen <- nchar(c4SamTable$read_seq)
    #c4SamTable[reference_name== 'hg38_knownGene_C4INSREG', ref_pos := ref_pos + 2860]
    
    if( nrow(c4SamTable) < 1000 ){
      next
    }
    
    ## Set the startPos and endPos for all reads
    c4SamTable[,startPos := mapply(run.setStartPos, locus, ref_pos)]
    c4SamTable[,endPos := mapply(run.setEndPos, locus, ref_pos, readLen)]
    
    ## Remove any positions that end up as NA (usually these are from reads that extend beyond the reference)
    c4SamTable <- c4SamTable[!is.na(c4SamTable$endPos)]
    c4SamTable <- c4SamTable[!is.na(c4SamTable$startPos)]
    
    ## Convert the SAM flags for all of the C4 aligned reads into an easily interpretable table
    bigSamFlagTable <- sapply(c4SamTable$sam_flag,samtable.flag_convert)
    
    c4SamTable <- c4SamTable[!unlist(bigSamFlagTable['notPrimaryAlignment',])]
    
    cat('\n\t\tConverting C4 aligned reads to read tables')
    c4SamTable[, 'readTable' := alleleSetup.read_formatter( read_name, read_seq, ref_pos, locus, cigar_string), by=seq_len(nrow(c4SamTable)) ]
    cat('\n\t\tDONE')
    
    c4BuildDF <- c4.samDT_to_depthDF(c4SamTable, resultsDirectory, currentSample$name, c4AlleleDF)
    
    ins.coordVect <- which(c4BuildDF['INS',] > 0)
    del.coordVect <- setdiff( which(c4BuildDF['.',] > 0), inverseDeletionIndexList$C4B )
    
    #importantRegions.cols <- c(14260, 14263, 14271, 14274, 14276, 14716, 14721, 14730, 14731, 15023, 15024)
    ## Phasing
    #nonINS.coord.vect <- setdiff( colnames(c4BuildDF), 2862:9226 )
    #hetCols <- which( apply( c4BuildDF >= minDP, 2, function(x) sum(x) > 1 ) )
    #if( length(hetCols) > 1 ){
    #  hetPos.toPhase.vect <- intersect(hetCols, nonINS.coord.vect)
    #  
    #  if( length(hetPos.toPhase.vect) > 1){
    #    hetSnpDF <- as.data.frame( apply( c4BuildDF[,hetPos.toPhase.vect], 2, function(x) which(x >= minDP)[1:2] ), check.names=F )
    #    hetPos.phased.list <- c4.snp_phaser(firstC4SamTable, secondC4SamTable, hetSnpDF, nucConv=T)
    #    c4.record_phasing(hetPos.phased.list, currentSample$name, phase.path)
    #  }
    #}
    
    #hetPos.phasedList <- c4.snp_phaser_2(firstC4SamTable, secondC4SamTable, importantRegions.cols)
    #c4.record_phasing(hetPos.phasedList, currentSample$name, phase.path)
    
    cat('\nPlotting C4 alignments.')
    #### C4 plotting
    ## Initialize color scheme for plots
    pal1 <- c("black","purple","green","blue","orange","brown","red")
    pal1 <- setNames(pal1, c('total','ANuc','TNuc','GNuc','CNuc','delNuc','insNuc'))
    
    #c4BuildList[[currentSample$name]] <- c4BuildDF
    
    miniDF <- c4BuildDF
    
    depthVect <- unlist(apply(miniDF,2,sum))
    
    p1 <- plot_ly(colors=pal1) %>%
      add_trace(x=as.integer(1:ncol(miniDF)),
                y=depthVect,
                mode='markers',type='scatter',name='total',color='total') %>%
      add_trace(x=as.integer(1:ncol(miniDF)),
                y=unlist(miniDF[1,]),
                mode='markers',type='scatter',name='A_nuc',color='ANuc') %>%
      add_trace(x=as.integer(1:ncol(miniDF)),
                y=unlist(miniDF[2,]),
                mode='markers',type='scatter',name='T_nuc',color='TNuc') %>%
      add_trace(x=as.integer(1:ncol(miniDF)),
                y=unlist(miniDF[3,]),
                mode='markers',type='scatter',name='C_nuc',color='CNuc') %>%
      add_trace(x=as.integer(1:ncol(miniDF)),
                y=unlist(miniDF[4,]),
                mode='markers',type='scatter',name='G_nuc',color='GNuc') %>%
      add_trace(x=as.integer(1:ncol(miniDF)),
                y=unlist(miniDF[5,]),
                mode='markers',type='scatter',name='del_nuc',color='delNuc') %>%
      add_trace(x=as.integer(1:ncol(miniDF)),
                y=unlist(miniDF[6,]),
                mode='markers',type='scatter',name='ins_nuc',color='insNuc')
    
    
    #print(p1)
    
    ## Save each plot
    htmlwidgets::saveWidget(p1, file=file.path(resultsDirectory,paste0(currentSample$name,'_c4_plot.html')))
    
    #### Various depth calculations
    
    cat('\nCalculating various depths.')
    ### Calculating TNXB depth
    ## Isolate TNXB aligned reads as a data table
    tnxbSamTable <- samTable[locus == 'TNXB']
    ## Convert the SAM flags to table format
    tnxbSamFlagTable <- sapply(tnxbSamTable$sam_flag, samtable.flag_convert)
    
    ## Isolate first-in-pair and second-in-pair reads as their own tables
    firstTnxbSamTable <- tnxbSamTable[unlist(tnxbSamFlagTable['firstInPair',])]
    secondTnxbSamTable <- tnxbSamTable[unlist(tnxbSamFlagTable['secondInPair',])]
    
    ## Calculate total TNXB depth as a function of the number of aligned reads multiplied by the median of the read lengths
    firstBaseCovNum <- length(unique(firstTnxbSamTable$read_name)) * median(as.integer(nchar(firstTnxbSamTable$read_seq)))
    secondBaseCovNum <- length(unique(secondTnxbSamTable$read_name)) * median(as.integer(nchar(secondTnxbSamTable$read_seq)))
    
    ## Calculate TNXB average depth by adding up both total first-in-pair and second-in-pair depth, then dividing by the number of bases in the TNXB reference
    medianTnxbDepthNum <- (firstBaseCovNum + secondBaseCovNum) / 68451
    
    ### Calculating C4 depth
    ## Caluculate average overall C4 depth over a stable region (positions 12k-17k in the long reference)
    medianc4DepthNum <- median(unlist(apply(c4BuildDF[,c(12000:17000)], 2, sum)))
    medianc4DepthNum <- median(unlist(apply(c4BuildDF[,c(17500:20000)], 2, sum)))    
    
    c4aMean <- mean(c(c4BuildDF[3,'E26_129'],c4BuildDF[4,'E26_132'],c4BuildDF[2,'E26_140'],c4BuildDF[4,'E26_143'],
                      c4BuildDF[3,'E26_145'], c4BuildDF[4,'E28_111'],c4BuildDF[2,'E28_116'],c4BuildDF[2,'E28_125'],c4BuildDF[3,'E28_126']))
    
    c4bMean <- mean(c(c4BuildDF[2,'E26_129'],c4BuildDF[3,'E26_132'],c4BuildDF[1,'E26_140'],c4BuildDF[3,'E26_143'],
                      c4BuildDF[2,'E26_145'], c4BuildDF[3,'E28_111'],c4BuildDF[3,'E28_116'],c4BuildDF[4,'E28_125'],c4BuildDF[4,'E28_126']))
    
    c4aSnpGroup1Mean <- mean(c(c4BuildDF[3,'E26_129'],c4BuildDF[4,'E26_132'],c4BuildDF[2,'E26_140'],c4BuildDF[4,'E26_143'],
                               c4BuildDF[3,'E26_145']))
    
    c4aSnpGroup2Mean <- mean(c(c4BuildDF[4,'E28_111'],c4BuildDF[2,'E28_116'],c4BuildDF[2,'E28_125'],c4BuildDF[3,'E28_126']))
    
    c4bSnpGroup1Mean <- mean(c(c4BuildDF[2,'E26_129'],c4BuildDF[3,'E26_132'],c4BuildDF[1,'E26_140'],c4BuildDF[3,'E26_143'],
                               c4BuildDF[2,'E26_145']))
    
    c4bSnpGroup2Mean <- mean(c(c4BuildDF[3,'E28_111'],c4BuildDF[3,'E28_116'],c4BuildDF[4,'E28_125'],c4BuildDF[4,'E28_126']))
    
    c4aLongSnp <- c4BuildDF[4,'I9_4288']
    
    c4bLongSnp <- c4BuildDF[1,'I9_4288']
    
    c4E29INS <- c4BuildDF[6,'E29_16']
    c4E29NORM <- c4BuildDF[4,'E29_16']
    
    #calculating median depth of deletion positions over the insertion region
    deletionDepth <- median(c(as.numeric(c4BuildDF[5,2861:9227])))
    
    #calculating median depth of the insertion region
    nonDeletionDepth <- mean( apply( c4BuildDF[1:4,2861:9227], 2, sum) )
    
    wgsC4Copy <- (medianc4DepthNum/TNXBMedian)*2
    resultsDF[currentSample$name,c('median_mhc_c4', 'median_mhc_tnxb', 'wgs_c4_copy')] <- c(c4Median, TNXBMedian, wgsC4Copy)
    
    #resultsDF[currentSample$name,'c4_exon_del'] <- paste0( c4.delPosVect, collapse='.' )
    resultsDF[currentSample$name,c('c4_exon29_ins','c4_exon29_norm')] <- c(c4E29INS,c4E29NORM)
    
    resultsDF[currentSample$name,c('median_c4mid','median_c4del','mean_c4ins','median_tnxb','mean_c4a_g1', 'mean_c4a_g2','mean_c4b_g1', 'mean_c4b_g2',
                                   'mean_c4a', 'mean_c4b', 'c4aL', 'c4bL')] <- c(medianc4DepthNum, deletionDepth, nonDeletionDepth, medianTnxbDepthNum,
                                                                                 c4aSnpGroup1Mean, c4aSnpGroup2Mean, c4bSnpGroup1Mean, c4bSnpGroup2Mean, c4aMean, c4bMean, c4aLongSnp, c4bLongSnp)
    
    #write.csv(c4BuildDF, file=file.path(resultsDirectory,paste0(projectName,'_c4_dp.csv')))                                  
    write.csv(resultsDF, file=file.path(resultsDirectory,paste0(projectName,'_c4_depth.csv')))
    file.remove(currentSample[['mhcBamPath']])
    #file.remove(currentSample[['bamPath']])
    file.remove(currentSample[['sortedBamPath']])
    file.remove(currentSample[['allDepthPath']])
  }
  

  return(NULL)
}
#write.csv(resultsDF, file=file.path(resultsDirectory,paste0(projectName,'_c4_depth.csv')))

c4.run_script( sampleList, projectName, outDir$path, threads, maxReadThreshold, minDP )

out_d <- normalizePath(outDir$path)
setwd(owd)
tar(output_f,out_d , compression = 'gzip', tar="tar")
unlink(out_d,recursive=TRUE)
