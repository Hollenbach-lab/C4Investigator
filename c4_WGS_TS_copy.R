library(data.table)
library(methods)
library(stringr)
library(plotly)

setwd('/home/LAB_PROJECTS/scripts')

source('function_base/general_functions.R')
source('function_base/c4Copy_functions.R')
source('function_base/generalAlignment_functions.R')

sampleDirectory  <- '/home/LAB_PROJECTS/MHC_PROJECT/1_raw_data/US_EURO_descendant/'
projectName <- 'C4_PROJECT'
resultsDirectory <- 'new_run_test'
fastqPattern <- 'fastq'
threads <- 40
maxReadThreshold <- 50000
minDP <- 6

c4.run_script <- function( sampleDirectory, projectName, resultsDirectory, fastqPattern, threads, maxReadThreshold, minDP ){

  sampleDirectory <- normalizePath(sampleDirectory, mustWork=T)
  resourcesDirectory <- normalizePath('/home/LAB_PROJECTS/reference/C4_copy_reference/', mustWork=T)
  mhcResourcesDirectory <- normalizePath('/home/LAB_PROJECTS/reference/bowtie2_index/chr6/', mustWork=T)
  
  baseOutputDir <- file.path('/home/LAB_PROJECTS',projectName,'output/')
  
  ## Setting up the results directory
  resultsDirectory <- file.path(baseOutputDir,resultsDirectory)
  dir.create(resultsDirectory) ## This function will output an unconsequential warning if the directory already exists
  resultsDirectory <- normalizePath(resultsDirectory, mustWork=T) ## Make sure the results directory exists
  
  ## Setting up C4 reference resources
  referencePath <- file.path(resourcesDirectory,'all_onelines_oneDel_bShort') ## C4 alignment reference
  mhcReferencePath <- file.path(mhcResourcesDirectory,'mhc')
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
  
  ## Read in probe CSV
  c4.probeDF <- read.csv(file=file.path('/home/LAB_PROJECTS/scripts/C4_exon_frame_shift_mut_probes.csv'),check.names=F,
                         stringsAsFactors = F)
  row.names(c4.probeDF) <- c4.probeDF$ProbeID
  
  ## Build up deletion index lists, which are needed for converting read coordinates between Along and Bshort
  inverseDeletionIndexList <- build.c4_inverse_deletion_index_list(c4AlleleDF)
  deletionIndexList <- build.c4_deletion_index_list(c4AlleleDF)
  
  #### These sets of functions use the deletion index lists to convert read alignment coordinates to a common coordinate
  ## This function sets the universal read start position
  run.setStartPos <- function(ref_name, ref_pos){
    return(deletionIndexList[[ref_name]][ref_pos])
  }
  ## This function sets the universal read end position
  run.setEndPos <- function(ref_name, ref_pos, currentReadLen){
    return(deletionIndexList[[ref_name]][ref_pos+currentReadLen-1])
  }
  
  ## Read in fastq samples from the sampleDirectory
  sampleList <- sequence.paired_sample_objects(sampleDirectory, fastqPattern, resultsDirectory)
  
  #sampleList <- c(sampleList1, sampleList2, sampleList3, sampleList4)
  
  #'specific to EUR FTD data'
  #sampleList <- c(sampleList1, sampleList2, sampleList3)
  #ftdIDFrame <- read.table('/home/LAB_PROJECTS/UCSF_MAC/EUR_CTL_leftover_after_lastBatch_sampleIDs.txt',stringsAsFactors = F)
  #ftdIDVect <- ftdIDFrame$V1
  
  #toRunList <- list()
  #for(sample in sampleList){
  #  sampleID <- tstrsplit(sample$name, '_', fixed=T)[[1]]
  #  
  #  if(sampleID %in% ftdIDVect){
  #    toRunList[[sample$name]] <- sample
  #  }
  #}
  
  #sampleList <- toRunList
  
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
  
  indel.path <- c4.initialize_indel_file(resultsDirectory)
  
  phase.path <- c4.initialize_phasing_file(resultsDirectory)
  
  ## Initialize a list for storing C4 alignments
  c4BuildList <- list()
  
  ## Initialize a dataframe for storing various depth values
  resultsDF <- data.frame(matrix(0,nrow=length(sampleList),ncol=18),row.names = names(sampleList),stringsAsFactors=F)
  colnames(resultsDF) <- c('median_c4mid','median_c4del','mean_c4ins','median_tnxb','mean_c4a_g1', 'mean_c4a_g2','mean_c4b_g1', 'mean_c4b_g2',
                           'mean_c4a', 'mean_c4b', 'c4aL', 'c4bL', 'median_mhc_c4', 'median_mhc_tnxb', 'mhc_c4_copy','c4_exon_del','c4_exon29_ins','c4_exon29_norm')
  
  cffCountDF <- data.frame(matrix(0,nrow=length(sampleList),ncol=nrow(c4.probeDF)), row.names=names(sampleList), stringsAsFactors = F)
  colnames(cffCountDF) <- rownames( c4.probeDF )
  
  for(currentSample in sampleList){
    
    ## If there is not a BAM file for this sample, then perform an alignment and convert the output to a BAM file
    if(!file.exists(currentSample$bamPath)){
      currentSample <- bowtie2.c4_alignment(bowtie2, referencePath, threads, currentSample, resultsDirectory)
      
      currentSample <- samtools.sam_to_bam(samtools, currentSample, resultsDirectory, threads)
      
      currentSample <- bowtie2.mhc.c4_alignment(bowtie2, mhcReferencePath, threads, currentSample, resultsDirectory)
      
      currentSample <- samtools.mhc.sam_to_bam(samtools, currentSample, resultsDirectory, threads)
    }else{
      ## If there is a BAM file for this sample, then convert the BAM file to a SAM file so it can be read in
      currentSample <- samtools.bam_to_sam(samtools, currentSample, resultsDirectory, threads)
      
      currentSample <- samtools.mhc.bam_to_sam(samtools, currentSample, resultsDirectory, threads)
    }
  
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
    samTable <- samfile.read_whole_genome_sam(currentSample$samPath, headerLineCountInt)
    
    file.remove(currentSample$samPath) ## Remove the SAM file to save space
    file.remove(currentSample$mhcSamPath) ## Remove the SAM file to save space
    
    ## Count CFF matches
    c4.mutCountList <- c4.count_cff_matches(currentSample, c4.probeDF, samTable)
    cffCountDF[currentSample$name,names(c4.mutCountList)] <- c4.mutCountList
    write.csv(cffCountDF, file=file.path(resultsDirectory,paste0(projectName,'_c4_cffCount.csv')))
  
    ## Convert the alignment scores into integers
    alignmentScoreList <- as.integer(tstrsplit(samTable$'12', ':', fixed=TRUE)[[3]])
  
    ## Save alignment scores to their own columns
    cat('\nSetting alignment scores.')
    samTable[,alignment_score := alignmentScoreList]
  
    samTable$reference_name <- tstrsplit(samTable$reference_name, '_')[[3]]
  
    ## Use the reference name conversion list to apply common gene names
    samTable$locus <- unlist(lapply(samTable$reference_name, function(x){referenceKeyList[[x]]}))
  
    ## Pull out an index of reads that aligned to C4
    C4ReadsPosVect <- grep('C4',samTable$locus,fixed=T)
  
    ## Initialize a column of the data table for storing whether the read aligned to C4 or not
    samTable$isC4 <- F
    samTable$isC4[C4ReadsPosVect] <- T ## Set the C4 bool to T for all C4 indexed reads
  
    ## Pull out the subset of C4 aligned reads as its own data table
    c4SamTable <- samTable[isC4 == T][locus %in% c('C4A','C4B')]
    
    if( nrow(c4SamTable) < 1000 ){
      next
    }
  
    ## Intialize and set a column for read lengths
    c4SamTable$readLen <- nchar(c4SamTable$read_seq)
  
    ## Set the startPos and endPos for all reads
    c4SamTable[,startPos := mapply(run.setStartPos, locus, ref_pos)]
    c4SamTable[,endPos := mapply(run.setEndPos, locus, ref_pos, readLen)]
  
    ## Remove any positions that end up as NA (usually these are from reads that extend beyond the reference)
    c4SamTable <- c4SamTable[!is.na(c4SamTable$endPos)]
    c4SamTable <- c4SamTable[!is.na(c4SamTable$startPos)]
  
    ## Convert the SAM flags for all of the C4 aligned reads into an easily interpretable table
    bigSamFlagTable <- sapply(c4SamTable$sam_flag,samtable.flag_convert)
  
    ## Isolate all first-in-pair reads (reads that aligned as the first half of a paired-end read pair)
    firstC4SamTable <- c4SamTable[unlist(bigSamFlagTable['firstInPair',])]
    ## Isolate all second-in-pair reads
    secondC4SamTable <- c4SamTable[unlist(bigSamFlagTable['secondInPair',])]
  
    ## Sort the paired-end sam tables by ascending read names and descending alignment scores
    firstC4SamTable <- firstC4SamTable[order(read_name, -alignment_score)]
    secondC4SamTable <- secondC4SamTable[order(read_name, -alignment_score)]
  
    ## Keep the alignments with the max score for each read
    secondC4SamTable <- secondC4SamTable[secondC4SamTable[,alignment_score == max(alignment_score),by=.(read_name)]$V1]
    firstC4SamTable <- firstC4SamTable[firstC4SamTable[,alignment_score == max(alignment_score),by=.(read_name)]$V1]
  
    ## Remove duplicate read alignments
    firstC4SamTable <- unique(firstC4SamTable, by=c('read_name'))
    secondC4SamTable <- unique(secondC4SamTable, by=c('read_name'))
  
    ## Remove reads with 'N' characters
    firstC4SamTable <- firstC4SamTable[grep('N',firstC4SamTable$read_seq,fixed=T,invert=T)]
    secondC4SamTable <- secondC4SamTable[grep('N',secondC4SamTable$read_seq,fixed=T,invert=T)]
  
    
    ## Identification of INSERTIONS
    first.exonINS.list <- c4.exon_ins_pos_list(firstC4SamTable, c4.exonCoordVect)
    second.exonINS.list <- c4.exon_ins_pos_list(secondC4SamTable, c4.exonCoordVect)
    total.exonINS.list <- c( first.exonINS.list, second.exonINS.list )
    c4.record_exon_ins(total.exonINS.list, resultsDirectory, currentSample$name, indel.path)
    
    
    
    currentDelIndex <<- inverseDeletionIndexList
    
    cat('\nFormatting paired C4 aligned reads.')
    ## Apply the read_formatter function to the overlapping reads (creates a named vector of the read, adds in deletion positions)
    firstC4SamTable[,read_table:=mapply(c4.read_formatter, read_name, read_seq, startPos, endPos, locus, cigar_string) ]
    secondC4SamTable[,read_table:=mapply(c4.read_formatter, read_name, read_seq, startPos, endPos, locus, cigar_string) ]
    
    ## Initialize a dataframe for storing the read information
    c4BuildDF <- build.c4_nuc_frame(c4AlleleDF)
    
    maxCol <- as.integer( colnames(c4BuildDF)[ncol(c4BuildDF)] )
    
    toRemove.index <- which( sapply(firstC4SamTable$read_table, function(x) max(as.integer(names(x))) > maxCol) )
    if(length(toRemove.index) > 0){
      firstC4SamTable <- firstC4SamTable[!toRemove.index,]
    }
    
    toRemove.index <- which( sapply(secondC4SamTable$read_table, function(x) max(as.integer(names(x))) > maxCol) )
    if(length(toRemove.index) > 0){
      secondC4SamTable <- secondC4SamTable[!toRemove.index,]
    }
  
    cat('\nBuilding C4 alignments.')
  
    for(currentTable in firstC4SamTable$read_table){
      ## Convert the nucleotide table into an integer table (nucListConv is in gc_functions.R)
      currentTableConverted <- lapply(currentTable, function(x) nucListConv[[x]])
  
      ## Fill in the build for the current locus using the converted table
      for(nucInt in unique(unlist(currentTableConverted))){
        coordsMatchingNucInt <- names(currentTableConverted[currentTableConverted == nucInt])
        c4BuildDF[as.character(nucInt),coordsMatchingNucInt] <- c4BuildDF[as.character(nucInt),coordsMatchingNucInt] + 1
      }
    }
  
    for(currentTable in secondC4SamTable$read_table){
      ## Convert the nucleotide table into an integer table (nucListConv is in gc_functions.R)
      currentTableConverted <- lapply(currentTable, function(x) nucListConv[[x]])
  
      ## Fill in the build for the current locus using the converted table
      for(nucInt in unique(unlist(currentTableConverted))){
        coordsMatchingNucInt <- names(currentTableConverted[currentTableConverted == nucInt])
        c4BuildDF[as.character(nucInt),coordsMatchingNucInt] <- c4BuildDF[as.character(nucInt),coordsMatchingNucInt] + 1
      }
    }
  
    ## DEL identification and recording
    c4ExonDelIndex <- which( c4BuildDF[5,c4.exonCoordVect] > minDP )
  
    c4.delPosVect <- c(0)
    if( length(c4ExonDelIndex) > 0 ){
      c4.delPosVect <- c4.exonCoordVect[c4ExonDelIndex]
    }
  
    exonDEL.index <- which( c4BuildDF[5,c4.exonCoordVect] > 0 )
    if( length(exonDEL.index) > 0 ){
      exonDELpos.vect <- c4.exonCoordVect[exonDEL.index]
      exonDELdp.vect <- c4BuildDF[5,exonDELpos.vect]
  
      lapply(1:length(exonDELpos.vect), function(x){
        delPos <- exonDELpos.vect[x]
        delDP <- exonDELdp.vect[x]
        textStr <- paste0(currentSample$name,'\t',delPos,'\t',delDP,'\tDEL\n')
        cat(textStr, file=indel.path, append=T)
      })
    }
    
    importantRegions.cols <- c(14260, 14263, 14271, 14274, 14276, 14716, 14721, 14730, 14731, 15023, 15024)
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
    
    hetPos.phasedList <- c4.snp_phaser_2(firstC4SamTable, secondC4SamTable, importantRegions.cols)
    c4.record_phasing(hetPos.phasedList, currentSample$name, phase.path)
    
    cat('\nPlotting C4 alignments.')
    #### C4 plotting
    ## Initialize color scheme for plots
    pal1 <- c("black","purple","green","blue","orange","brown")
    pal1 <- setNames(pal1, c('total','ANuc','TNuc','GNuc','CNuc','delNuc'))
  
    #c4BuildList[[currentSample$name]] <- c4BuildDF
  
    miniDF <- c4BuildDF
  
    depthVect <- unlist(apply(miniDF,2,sum))
  
    p1 <- plot_ly(colors=pal1) %>%
      add_trace(x=as.integer(colnames(miniDF)),
                y=depthVect,
                mode='markers',type='scatter',name='total',color='total') %>%
      add_trace(x=as.integer(colnames(miniDF)),
                y=unlist(miniDF[1,]),
                mode='markers',type='scatter',name='A_nuc',color='ANuc') %>%
      add_trace(x=as.integer(colnames(miniDF)),
                y=unlist(miniDF[2,]),
                mode='markers',type='scatter',name='T_nuc',color='TNuc') %>%
      add_trace(x=as.integer(colnames(miniDF)),
                y=unlist(miniDF[3,]),
                mode='markers',type='scatter',name='C_nuc',color='CNuc') %>%
      add_trace(x=as.integer(colnames(miniDF)),
                y=unlist(miniDF[4,]),
                mode='markers',type='scatter',name='G_nuc',color='GNuc') %>%
      add_trace(x=as.integer(colnames(miniDF)),
                y=unlist(miniDF[5,]),
                mode='markers',type='scatter',name='del_nuc',color='delNuc')
  
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
  
    c4aMean <- mean(c(c4BuildDF[3,'14260'],c4BuildDF[4,'14263'],c4BuildDF[2,'14271'],c4BuildDF[4,'14274'],
                      c4BuildDF[3,'14276'], c4BuildDF[4,'14716'],c4BuildDF[2,'14721'],c4BuildDF[2,'14730'],c4BuildDF[3,'14731']))
  
    c4bMean <- mean(c(c4BuildDF[2,'14260'],c4BuildDF[3,'14263'],c4BuildDF[1,'14271'],c4BuildDF[3,'14274'],
                      c4BuildDF[2,'14276'], c4BuildDF[3,'14716'],c4BuildDF[3,'14721'],c4BuildDF[4,'14730'],c4BuildDF[4,'14731']))
  
    c4aSnpGroup1Mean <- mean(c(c4BuildDF[3,'14260'],c4BuildDF[4,'14263'],c4BuildDF[2,'14271'],c4BuildDF[4,'14274'],
                               c4BuildDF[3,'14276']))
  
    c4aSnpGroup2Mean <- mean(c(c4BuildDF[4,'14716'],c4BuildDF[2,'14721'],c4BuildDF[2,'14730'],c4BuildDF[3,'14731']))
  
    c4bSnpGroup1Mean <- mean(c(c4BuildDF[2,'14260'],c4BuildDF[3,'14263'],c4BuildDF[1,'14271'],c4BuildDF[3,'14274'],
                               c4BuildDF[2,'14276']))
  
    c4bSnpGroup2Mean <- mean(c(c4BuildDF[3,'14716'],c4BuildDF[3,'14721'],c4BuildDF[4,'14730'],c4BuildDF[4,'14731']))
  
    c4aLongSnp <- c4BuildDF[4,'6873']
  
    c4bLongSnp <- c4BuildDF[1,'6873']
  
    c4E29INS <- c4BuildDF[4,'15028']
    c4E29NORM <- c4BuildDF[3,'15028']
  
    #calculating median depth of deletion positions over the insertion region
    deletionDepth <- median(c(as.numeric(c4BuildDF[5,as.character(2861:9227)])))
  
    #calculating median depth of the insertion region
    nonDeletionDepth <- mean( apply( c4BuildDF[1:4,as.character(2861:9227)], 2, sum) )
  
    wgsC4Copy <- (medianc4DepthNum/TNXBMedian)*2
    resultsDF[currentSample$name,c('median_mhc_c4', 'median_mhc_tnxb', 'wgs_c4_copy')] <- c(c4Median, TNXBMedian, wgsC4Copy)
  
    resultsDF[currentSample$name,'c4_exon_del'] <- paste0( c4.delPosVect, collapse='.' )
    resultsDF[currentSample$name,c('c4_exon29_ins','c4_exon29_norm')] <- c(c4E29INS,c4E29NORM)
  
    resultsDF[currentSample$name,c('median_c4mid','median_c4del','mean_c4ins','median_tnxb','mean_c4a_g1', 'mean_c4a_g2','mean_c4b_g1', 'mean_c4b_g2',
                                   'mean_c4a', 'mean_c4b', 'c4aL', 'c4bL')] <- c(medianc4DepthNum, deletionDepth, nonDeletionDepth, medianTnxbDepthNum,
                                        c4aSnpGroup1Mean, c4aSnpGroup2Mean, c4bSnpGroup1Mean, c4bSnpGroup2Mean, c4aMean, c4bMean, c4aLongSnp, c4bLongSnp)
  
    write.csv(resultsDF, file=file.path(resultsDirectory,paste0(projectName,'_c4_depth.csv')))
  }
  return(NULL)
}
#write.csv(resultsDF, file=file.path(resultsDirectory,paste0(projectName,'_c4_depth.csv')))



