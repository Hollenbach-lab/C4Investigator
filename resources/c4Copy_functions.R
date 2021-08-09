

## This function runs a bowtie2 alignment looking for exact matches to all KIR references in subsetKirReference
bowtie2.c4_alignment <- function(bowtie2_command, reference_index, threads, current_sample, resultsDirectory){
  
  ## Intitialize an output path for the SAM file
  current_sample[['samPath']] <- file.path(resultsDirectory,paste0(current_sample$name,'.sam'))
  
  ## Building up the run command
  optionsCommand <- c(paste0('-x ',reference_index),
                      '-5 0', '-3 6', '-N 0', '--end-to-end', paste0('-p ',threads), '--score-min "L,0,-0.10"',
                      '-I 75', '-X 1000',
                      paste0('-1 ',current_sample$fastq1path),
                      paste0('-2 ',current_sample$fastq2path),
                      '--no-unal','-a','--mp 2,2', '--rdg 1,1', '--rfg 1,1',
                      paste0('-S ',current_sample$samPath))
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2_command,optionsCommand)
  output.sampleAlign <- system2(bowtie2_command, optionsCommand, stdout=T, stderr=T)
  
  if(!is.null(attributes(output.sampleAlign))){
    cat('\nBowtie2 failed, retrying alignment...')
    output.sampleAlign <- system2(bowtie2_command, optionsCommand, stdout=T, stderr=T)
  }
  
  check.system2_output(output.sampleAlign, 'bowtie2 gc alignment failed')
  
  ## Print the bowtie2 output
  cat('\n',paste0(output.sampleAlign, collapse='\n'))
  
  ## Check to make sure the SAM file actually exists
  current_sample[['samPath']] <- normalizePath(current_sample$samPath, mustWork=T)
  
  cat('\n\nSuccessfully aligned',current_sample$name,'to',reference_index)
  
  return(current_sample)
}

## This function runs a bowtie2 alignment looking for exact matches to all KIR references in subsetKirReference
bowtie2.mhc.c4_alignment <- function(bowtie2_command, reference_index, threads, current_sample, resultsDirectory){
  
  ## Intitialize an output path for the SAM file
  current_sample[['mhcSamPath']] <- file.path(resultsDirectory,paste0(current_sample$name,'.mhc.sam'))
  
  ## Building up the run command
  optionsCommand <- c(paste0('-x ',reference_index),
                      '-5 0', '-3 6', '-N 0', '--end-to-end', paste0('-p ',threads), '--score-min "L,0,-0.10"',
                      '-I 75', '-X 1000',
                      paste0('-1 ',current_sample$fastq1path),
                      paste0('-2 ',current_sample$fastq2path),
                      '--no-unal','-a','--mp 2,2', '--rdg 1,1', '--rfg 1,1',
                      paste0('-S ',current_sample$mhcSamPath))
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2_command,optionsCommand)
  output.sampleAlign <- system2(bowtie2_command, optionsCommand, stdout=T, stderr=T)
  
  if(!is.null(attributes(output.sampleAlign))){
    cat('\nBowtie2 failed, retrying alignment...')
    output.sampleAlign <- system2(bowtie2_command, optionsCommand, stdout=T, stderr=T)
  }
  
  check.system2_output(output.sampleAlign, 'bowtie2 gc alignment failed')
  
  ## Print the bowtie2 output
  cat('\n',paste0(output.sampleAlign, collapse='\n'))
  
  ## Check to make sure the SAM file actually exists
  current_sample[['mhcSamPath']] <- normalizePath(current_sample$mhcSamPath, mustWork=T)
  
  cat('\n\nSuccessfully aligned',current_sample$name,'to',reference_index)
  
  return(current_sample)
}
## This function returns a list of dataframes of allele sequences found in the reference fasta
read.c4_dataframe_from_reference_fasta <- function(fasta_path, referenceKeyList){
  
  ## Make sure the fasta can be found
  fasta_path <- normalizePath(fasta_path, mustWork=T)
  
  ## Initialize a list to store the sequence strings from the file
  alleleSeqList <- list()
  
  ## Read in the fasta file and store the allele names and sequences
  for(currentLine in readLines(fasta_path)){
    alleleNameBool <- grepl('>',currentLine,fixed=T)
    
    if(alleleNameBool){
      alleleName <- strsplit(currentLine, '_',fixed=TRUE)[[1]][3]
      alleleName <- strsplit(alleleName,' ',fixed=T)[[1]][1]
      alleleName <- referenceKeyList[[alleleName]]
    }else{
      alleleSeq <- currentLine
      
      if(alleleName %in% c('C4A','C4B')){
        alleleSeqList[[alleleName]] <- alleleSeq
      }
    }
  }
  
  ## Initialize the output list
  output.alleleSeqDFList <- list()
  
  cat('\n\tProcessing C4 sequence...')
  
  largestAlleleSize <- max(sapply(alleleSeqList,nchar))
  minAlleleSize <- min(sapply(alleleSeqList,nchar))
  
  ## Make sure the largest and smallest allele are the same sizes
  if(largestAlleleSize != minAlleleSize){
    stop(fasta_path,' has alleles of different sizes. Cannot transform into dataframe.')
  }
  
  ## Initialize a dataframe for storing the allele sequence strings
  alleleSeqDF <- data.frame(matrix(0, length(alleleSeqList), largestAlleleSize),row.names=names(alleleSeqList),check.names=F,stringsAsFactors=F)
  
  ## For each allele for the current locus, input the sequence into the dataframe
  for(currentAlleleName in names(alleleSeqList)){
    alleleSeqDF[currentAlleleName,] <- strsplit(alleleSeqList[[currentAlleleName]],'')[[1]]
  }
  
  return(alleleSeqDF)
}

## This function builds a list of deletion indices for each kir allele
build.c4_inverse_deletion_index_list <- function(c4AlleleDF){
  cat('\nBuilding up a deletion index for faster lookup during read assignment.')
  
  ## Initialize a list for storing deletion position conversions
  deletionIndexList <- list()
  for(currentAllele in rownames(c4AlleleDF)){
    deletionIndexList[[currentAllele]] <- which(c4AlleleDF[currentAllele,] == '.')
  }
  return(deletionIndexList)
}

## This function builds a list of deletion indices for each kir allele
build.c4_deletion_index_list <- function(c4AlleleDF){
  cat('\nBuilding up a deletion index for faster lookup during read assignment.')
  
  ## Initialize a list for storing deletion position conversions
  deletionIndexList <- list()
  for(currentAllele in rownames(c4AlleleDF)){
    deletionIndexList[[currentAllele]] <- which(c4AlleleDF[currentAllele,] != '.')
  }
  return(deletionIndexList)
}

## This function initializes a list of dataframes for storing the assembled nucleotides for each kir locus
build.c4_nuc_frame <- function(c4AlleleDF){
  refSeqDF <- data.frame(matrix(0,5,ncol(c4AlleleDF)),check.names=F,stringsAsFactors=F)
  colnames(refSeqDF) <- as.character(colnames(c4AlleleDF))
  rownames(refSeqDF) <- c('1','2','3','4','5')
  return(refSeqDF)
}

## This function counts how many reads map to a unique locus or allele
c4.count_read_matches <- function(currentSample, samTable, alignedLocusVect, maxReadThreshold){
  
  ## Pull out the unique read names
  uniqueReadNames <- unique(samTable$read_name)
  
  ## Randomize the read name order
  set.seed(001) # just to make it reproducible
  randomUniqueReadNames <- sample(uniqueReadNames)
  
  ## Check if there are more reads than the threshold and take some out if so
  if(length(randomUniqueReadNames) > maxReadThreshold){
    randomUniqueReadNames <- randomUniqueReadNames[1:maxReadThreshold]
  }
  
  ## Initialize the list for storing reference matches
  uniqueLocusMatchList <- as.list(alignedLocusVect)
  names(uniqueLocusMatchList) <- alignedLocusVect
  uniqueLocusMatchList[alignedLocusVect] <- 0
  
  
  ## Initialize variables to check on the progress of the next for loop
  i = 1
  max_i = length(randomUniqueReadNames)
  checkAmount = ceiling(max_i/10)
  j=0
  
  ## Find the reference matches for each unique read name
  for(currentReadName in randomUniqueReadNames){
    
    ## Pull out the lines of the SAM file that match the current read name
    samSubsetTable <- samTable[read_name == currentReadName]
    
    if('C4ins' %in% samSubsetTable$locus){
      uniqueLocusMatchList['C4ins'] = uniqueLocusMatchList['C4ins'][[1]] + 1
      uniqueLocusMatchList['C4'] = uniqueLocusMatchList['C4'][[1]]+1
      
      ## Will display the percent completion every 10%
      i = i+1
      if(i%%checkAmount == 0){
        j = j+10
        cat(paste0(j,'% ', collapse = ''))
      }
      
      next
    }
    
    ## Find the best alignment score for this read
    maxAlignmentScore <- max(samSubsetTable$alignment_score)
    
    ## Pull out the current read name alignments that have the best alignment score
    samSubsetTable <- samSubsetTable[alignment_score == maxAlignmentScore]
    
    ## Pull out the unique locus names
    matchedLocusList <- samSubsetTable$locus
    
    ###### This section will count all locus matches (as opposed to only unique locus matches)
    #for(matchedLocus in matchedLocusList){
    #  uniqueLocusMatchList[matchedLocus] = uniqueLocusMatchList[matchedLocus][[1]] + 1
    #}
    #if('KIR2DL5A' %in% matchedLocusList & 'KIR2DL5B' %in% matchedLocusList){
    #  uniqueLocusMatchList['KIR2DL5'] = uniqueLocusMatchList['KIR2DL5'][[1]] + 1
    #}
    ###### /s
    
    ###### This section will count only unique locus matches (as opposed to all locus matches)
    ## If there is only 1 unique locus, then add 1 to the unique match count for that locus
    
    if(length(matchedLocusList) == 1){
      uniqueLocusMatchList[matchedLocusList] = uniqueLocusMatchList[matchedLocusList][[1]] + 1
      
      ## If there is only a single matching reference allele, then iterate the count of that allele
      #if(length(matchedAlleleList) == 1){
      #  uniqueAlleleMatchList[matchedAlleleList] = uniqueAlleleMatchList[matchedAlleleList][[1]] + 1
      #}
    }else if('C4A' %in% matchedLocusList & 'C4B' %in% matchedLocusList){
      uniqueLocusMatchList['C4'] = uniqueLocusMatchList['C4'][[1]] + 1
    }
    ###### /s
    
    ## Will display the percent completion every 10%
    i = i+1
    if(i%%checkAmount == 0){
      j = j+10
      cat(paste0(j,'% ', collapse = ''))
    }
  }
  cat('\n',i)
  cat("\n\nFinished counting!")
  
  ## Cutting the function short for now to test what a good threshold would be
  return(list(locusMatches = uniqueLocusMatchList))
  #locusRatio <- sapply(uniqueLocusMatchList, function(x) x/uniqueLocusMatchList$KIR3DL3)
}

## This functions counts how many reads match CFF probes
c4.count_cff_matches <- function(currentSample,probeDF,samTable){
  cat('\n\nCounting CFF probes for the current sample..')
  
  ## Initialize the list for cff probe match counts
  cffProbeMatchList <- as.list(probeDF$Name)
  names(cffProbeMatchList) <- probeDF$Name
  cffProbeMatchList[probeDF$Name] <- 0
  
  for(probeName in probeDF$Name){
    cffProbeMatchList[probeName] <- length(grep(probeDF[probeName,'Sequence'], samTable$read_seq, fixed=T)) + as.integer(cffProbeMatchList[[probeName]])
  }
  cat('\tFinished.')
  return(cffProbeMatchList)
}

## Format reads for snp phasing (fill in deletion positions with '.')
c4.read_formatter <- function(readName, readSeq, startPos, endPos, refAlleleName, cigarString){
  #cat('\n',readName,readSeq,startPos,endPos,refAlleleName, cigarString, length(currentDelIndex))
  cigarOperationVect <- c('M','I','D')
  
  cigarCheckVect <- c(cigarOperationVect,'1','2','3','4','5','6','7','8','9','0')
  
  ## Format the start position
  startPos <- as.numeric(startPos)
  
  ## Format the end position
  endPos <- as.numeric(endPos)
  
  ## Testing if there deletions found in the reference for the current read
  if(nchar(cigarString) > (nchar(nchar(readSeq)) + 1)){
    
    ## Split the cigar string up into a vector
    cigarVect <- strsplit(cigarString,'')[[1]]
    
    ## Make sure there are no weird codes that we havent accounted for
    if(!all(cigarVect %in% cigarCheckVect)){
      stop('c4.read_formatter found cigar operators that are not coded for ',readName)
    }
    
    ## Determine which position of the vector denote cigar operators
    operationPosVect <- which(cigarVect %in% cigarOperationVect)
    
    ## Initialize a list for storing the different operations and their positions
    operationList <- list()
    prevPos <- 1
    countInt <- 1
    
    ## Convert the cigar vect into the operation list
    for(operationPos in operationPosVect){
      operationList[[countInt]] <- cigarVect[prevPos:operationPos]
      prevPos <- operationPos+1
      countInt <- countInt+1
    }
    
    ## Initialize a start position for CIGAR operations
    prevPos <- 1
    subStrVect <- c()
    for(operationVect in operationList){
      ## Pull out the end position of the current operation
      operationEndPos <- as.integer(paste0(operationVect[1:(length(operationVect)-1)], collapse='')) + prevPos - 1
      
      ## Pull out the operation type
      operationTypeChar <- operationVect[length(operationVect)]
      
      ## If the end position is not an integer something went horribly wrong
      if(is.na(operationEndPos) | !(operationTypeChar %in% cigarOperationVect)){
        stop('c4.read_formatter something weird happened for ',readName)
      }
      
      if(operationTypeChar == 'M'){
        ## Add the good matches to the read string vect
        subStrVect <- c(subStrVect, substr(readSeq,prevPos,operationEndPos))
        prevPos <- operationEndPos + 1
      }else if(operationTypeChar == 'I'){
        
        ## If the insertion is in this specific position we want to keep it
        if((startPos + prevPos) == 14796 & refAlleleName == 'C4B'){
          subStrVect <- c(subStrVect, substr(readSeq,prevPos,operationEndPos))
        }
        
        ## Condition to catch TC frame-shift insertion
        #if((startPos + prevPos) >= 15023 & (startPos + prevPos) <= 15025 ){
        #  subStrVect <-  c(subStrVect, substr(readSeq,prevPos,operationEndPos))
        #}
        
        ## Building up all insertions now
        #if((startPos+prevPos) %in% c4.exonCoordVect | operationEndPos %in% c4.exonCoordVect){
        subStrVect <-  c(subStrVect, substr(readSeq,prevPos,operationEndPos))
        #}
        
        ## Otherwise we ignore it
        prevPos <- operationEndPos + 1
      }else if(operationTypeChar == 'D'){
        ## Create a string of '.' of the current operation length
        delStr <- paste0(replicate((operationEndPos - prevPos) + 1, '.'), collapse='')
        
        ## Add the del string to the read string vect
        subStrVect <- c(subStrVect, delStr)
      }
    }
    
    ## Collapse the subStrVect into the read sequence
    readSeq <- paste0(subStrVect,collapse='')
  }else if((endPos-startPos+1) > nchar(readSeq)){
    
   # cat('\n',refAlleleName, length(currentDelIndex))
    delIndexList <- which(startPos:endPos %in% currentDelIndex[[refAlleleName]])
    
    ## If the index list is empty, something went wrong
    if(length(delIndexList) == 0){
      stop('Something weird happened with the deletion index.')
    }
    
    ## Split apart the read sequence and add in deletion symbols
    subStringList <- c()
    for(delIndex in delIndexList){
      preString <- substr(readSeq, 1, delIndex-1)
      postString <- substr(readSeq, delIndex, nchar(readSeq))
      
      readSeq <- paste0(preString, '.', postString)
    }
  }
  
  ## Create a full index that cooresponds to the sequence nucleotides and where they go
  fullIndex <- startPos:(nchar(readSeq)+startPos-1)
  
  ## Turn the sequence string into a list
  seqList <- strsplit(readSeq,'')[[1]]
  
  if(length(seqList) != length(fullIndex)){
    cat('\n',startPos,endPos,refAlleleName)
  }
  
  names(seqList) <- fullIndex
  
  seqListList <- list(seqList)
  
  return(seqListList)
}

## This function runs a samtools sam to bam conversion
samtools.mhc.sam_to_bam <- function(samtools_command, currentSample, resultsDirectory, threads){
  
  ## Initialize an output path for the BAM file
  currentSample[['mhcBamPath']] <- file.path(resultsDirectory,paste0(currentSample$name,'.mhc.bam'))
  
  ## Building up the run command
  optionsCommand <- c('view',paste0('-@', threads),
                      currentSample$mhcSamPath, '-o', currentSample$mhcBamPath)
  
  cat('\n\n',samtools_command, optionsCommand)
  output.bamConv <- system2(samtools_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.bamConv, 'samtools sam to bam conversion failed')
  
  ## Print the conversion output
  cat('\n',paste0(output.bamConv), collapse='\n')
  
  ## Check to make sure the BAM file actually exists
  currentSample[['mhcBamPath']] <- normalizePath(currentSample$mhcBamPath, mustWork=T)
  
  cat('\n\nSuccessfully converted',currentSample$mhcSamPath,'to',currentSample$mhcBamPath)
  
  return(currentSample)
}

## This function runs a samtools bam to sam conversion
samtools.mhc.bam_to_sam <- function(samtools_command, currentSample, resultsDirectory, threads){
  
  ## Initialize an output path for the BAM file
  currentSample[['mhcBamPath']] <- file.path(resultsDirectory,paste0(currentSample$name,'.mhc.bam'))
  
  ## Intitialize an output path for the SAM file
  currentSample[['mhcSamPath']] <- file.path(resultsDirectory,paste0(currentSample$name,'.mhc.sam'))

  
  ## Building up the run command
  optionsCommand <- c('view',paste0('-@', threads),'-h',
                      currentSample[['mhcBamPath']], '-o', currentSample[['mhcSamPath']])
  
  cat('\n\n',samtools_command, optionsCommand)
  output.samConv <- system2(samtools_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.samConv, 'samtools bam to sam conversion failed')
  
  ## Print the conversion output
  cat('\n',paste0(output.samConv), collapse='\n')
  
  ## Check to make sure the BAM file actually exists
  currentSample[['mhcSamPath']] <- normalizePath(currentSample[['mhcSamPath']], mustWork=T)
  
  cat('\n\nSuccessfully converted', 
      currentSample[['mhcBamPath']], 'to', 
      currentSample[['mhcSamPath']] )
  
  return(currentSample)
}

## This function sorts a bam file
samtools.mhc.sort <- function(samtools_command, currentSample, resultsDirectory, threads){
  ## Initialize an output path for the BAM file
  currentSample[['sortedBamPath']] <- file.path(resultsDirectory,paste0(currentSample$name,'.sorted.bam'))
  
  ## Building up the run command
  optionsCommand <- c('sort',paste0('-@', threads),
                      currentSample$mhcBamPath, '-o', currentSample$sortedBamPath)
  
  cat('\n\n',samtools_command, optionsCommand)
  output.bamSort <- system2(samtools_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.bamSort, 'samtools BAM sorting failed')
  
  ## Print the conversion output
  cat('\n',paste0(output.bamSort), collapse='\n')
  
  ## Check to make sure the sorted BAM file exists
  currentSample[['sortedBamPath']] <- normalizePath(currentSample$sortedBamPath, mustWork=T)
  
  cat('\n\nSuccessfully sorted',currentSample$sortedBamPath)
  
  return(currentSample)
}



### This function attempts to phase het SNPs using unpaired and paired-end reads
c4.snp_phaser <- function(c4.1.samTable, c4.2.samTable, hetSnpDF, nucConv = T){
  
  cat('\nPhasing',ncol(hetSnpDF),'het SNPs')
  ## Initialize a list for storing the phased snps (new elements are added when phasing cannot be completed)
  phasedList <- list()
  
  ## Initialze a count for keeping track of what list element we are on
  i <- 1
  
  ## Initialize the previous snp position
  prevPos <- 0
  
  ## Loop over all het snps
  for(snpPos in colnames(hetSnpDF)){
    
    ## Pull out the nucleotides for the current position
    if(nucConv){
      snpVect <- names(nucListConv)[as.numeric(hetSnpDF[,snpPos])]
    }else{
      snpVect <- as.vector(currentLocusHetSnpDF[,snpPos])
    }
    
    ## If this is the first iteration, then initialize a dataframe for storing the phased Snps
    if(prevPos == 0){
      ## Data frame with two rows for storing strand1 and strand2 SNPs (s1 and s2 from different list elements do not match up)
      #phasedList[[i]] <- data.frame(matrix('',nrow=2,ncol=1),row.names=c('s1','s2'))
      phasedList[[i]] <- NA
      ## Name the dataframe colums based on the current snp position
      #colnames(phasedList[[i]]) <- snpPos
      
      ## Save the dataframe to the current list element
      #phasedList[[i]][,snpPos] <- snpVect
      prev.snp1 <- snpVect[1]
      prev.snp2 <- snpVect[2]
    }else{
      
      startingPairedReads <- c4.1.samTable[ startPos <= as.numeric(prevPos) ]$read_name
      startingPairedReads <- unique( c( startingPairedReads, c4.2.samTable[ startPos <= as.numeric(prevPos) ]$read_name ) )
      
      endingPairedReads <- c4.1.samTable[ endPos >= as.numeric(snpPos) ]$read_name
      endingPairedReads <- unique( c( endingPairedReads, c4.2.samTable[ endPos >= as.numeric(snpPos) ]$read_name ) )
      
      spanningPairedReads <- intersect(startingPairedReads, endingPairedReads)
      
      if(length(spanningPairedReads) == 0){
        #i <- i + 1
        #phasedList[[i]] <- data.frame(matrix('',nrow=2,ncol=1),row.names=c('s1','s2'))
        #colnames(phasedList[[i]]) <- snpPos
        #phasedList[[i]][,snpPos] <- snpVect
      }else{
        
        ## Pull out the prevPos snps
        #s1Snp1 <- phasedList[[i]]['s1',prevPos]
        #s2Snp1 <- phasedList[[i]]['s2',prevPos]
        s1Snp1 <- prev.snp1
        s2Snp1 <- prev.snp2
        
        ## Set the snpPos snps
        s1Snp2 <- snpVect[1]
        s2Snp2 <- snpVect[2]
        
        span.1.samTable <- unique( c4.1.samTable[ read_name %in% spanningPairedReads ], by=c('read_seq', 'read_name') )
        span.2.samTable <- unique( c4.2.samTable[ read_name %in% spanningPairedReads ], by=c('read_seq', 'read_name') )
        
        unique.spanningPairedReads <- unique(span.1.samTable$read_name, span.2.samTable$read_name)
        
        test.list <- list()
        test.iter <- 1
        
        for( currentRead in unique.spanningPairedReads){

          if( currentRead %in% span.1.samTable$read_name ){
            snp1.1 <- span.1.samTable[ read_name == currentRead ]$read_table[[1]][prevPos]
            snp2.1 <- span.1.samTable[ read_name == currentRead ]$read_table[[1]][snpPos]
          }else{
            snp1.1 <- NA
            snp2.1 <- NA
          }
          
          if( currentRead %in% span.2.samTable$read_name ){
            snp1.2 <- span.2.samTable[ read_name == currentRead ]$read_table[[1]][prevPos]
            snp2.2 <- span.2.samTable[ read_name == currentRead ]$read_table[[1]][snpPos]
          }else{
            snp1.2 <- NA
            snp2.2 <- NA
          }
          
          phasedSnps.vect <- c(snp1.1, snp2.1, snp1.2, snp2.2)
          phasedSnps.vect <- phasedSnps.vect[!is.na(phasedSnps.vect)]
          
          phasedSnps.vect <- phasedSnps.vect[ sort( unique( names(phasedSnps.vect) ) )]
          
          temp.names <- names( phasedSnps.vect )
          temp.output <- paste0( c(temp.names, phasedSnps.vect[temp.names]), collapse='^' )
          
          if(length( phasedSnps.vect) > 1){
            test.list[[test.iter]] <- temp.output
            test.iter <- test.iter+1
          }
        }
        
        phasedList[[i]] <- table( unlist( test.list ) )
        i <- i + 1
      }
      
      prev.snp1 <- snpVect[1]
      prev.snp2 <- snpVect[2]
    }
    
    ## Set the prevPos to be the current snpPos, then move on to the next snpPos
    prevPos <- snpPos
  }
  cat('\nDone.')
  return(phasedList)
}

### This function attempts to phase het SNPs using unpaired and paired-end reads
c4.snp_phaser_2 <- function(c4.1.samTable, c4.2.samTable, hetSnpVect){
  
  cat('\nPhasing important SNPs')
  ## Initialize a list for storing the phased snps (new elements are added when phasing cannot be completed)
  phasedList <- list()
  
  ## Initialze a count for keeping track of what list element we are on
  i <- 1
  
  ## Initialize the previous snp position
  prevPos <- 0
  
  ## Loop over all het snps
  for(snpPos in hetSnpVect){
    
    curPos <- as.character(snpPos)
    
    if(prevPos == 0){
      prevPos <- curPos
      next
    }
    
    startingPairedReads <- c4.1.samTable[ startPos <= as.numeric(prevPos) ]$read_name
    startingPairedReads <- unique( c( startingPairedReads, c4.2.samTable[ startPos <= as.numeric(prevPos) ]$read_name ) )
    
    endingPairedReads <- c4.1.samTable[ endPos >= as.numeric(curPos) ]$read_name
    endingPairedReads <- unique( c( endingPairedReads, c4.2.samTable[ endPos >= as.numeric(curPos) ]$read_name ) )
    
    spanningPairedReads <- intersect(startingPairedReads, endingPairedReads)
    
    if(length(spanningPairedReads) != 0){
      span.1.samTable <- unique( c4.1.samTable[ read_name %in% spanningPairedReads ], by=c('read_seq', 'read_name') )
      span.2.samTable <- unique( c4.2.samTable[ read_name %in% spanningPairedReads ], by=c('read_seq', 'read_name') )
      
      unique.spanningPairedReads <- unique(span.1.samTable$read_name, span.2.samTable$read_name)
      

      ## Isolate SNP calls from reads 
      pos1.read1.list <- sapply(unique.spanningPairedReads, function(x){
        result = tryCatch({
          as.character( span.1.samTable[ read_name == x ][,read_table[[1]][prevPos]] )
        }, error = function(e) {
          NA
        })
      })
      pos1.read1.list <- pos1.read1.list[!is.na( pos1.read1.list )]
      
      pos2.read1.list <- sapply(unique.spanningPairedReads, function(x){
        result = tryCatch({
          as.character( span.1.samTable[ read_name == x ][,read_table[[1]][curPos]] )
        }, error = function(e) {
          NA
        })
      })
      pos2.read1.list <- pos2.read1.list[!is.na( pos2.read1.list )]
      
      pos1.read2.list <- sapply(unique.spanningPairedReads, function(x){
        result = tryCatch({
          as.character( span.2.samTable[ read_name == x ][,read_table[[1]][prevPos]] )
        }, error = function(e) {
          NA
        })
      })
      pos1.read2.list <- pos1.read2.list[!is.na( pos1.read2.list )]
      
      pos2.read2.list <- sapply(unique.spanningPairedReads, function(x){
        result = tryCatch({
          as.character( span.2.samTable[ read_name == x ][,read_table[[1]][curPos]] )
        }, error = function(e) {
          NA
        })
      })
      pos2.read2.list <- pos2.read2.list[!is.na( pos2.read2.list )]
      
      
      ## Format phased calls into list
      matchedSnp.list <- sapply(unique.spanningPairedReads, function(x, prevPos, curPos){
        
        result1 <- tryCatch({
          paste(prevPos,curPos,pos1.read1.list[[x]],pos2.read1.list[[x]],sep='^')
        }, error= function(e){
          NA
        })
        
        result2 <- tryCatch({
          paste(prevPos,curPos,pos1.read1.list[[x]],pos2.read2.list[[x]],sep='^')
        }, error= function(e){
          NA
        })
        
        result3 <- tryCatch({
          paste(prevPos,curPos,pos1.read2.list[[x]],pos2.read2.list[[x]],sep='^')
        }, error= function(e){
          NA
        })
        
        out.vect <- c(result1, result2, result3)
        out.vect <- out.vect[!is.na(out.vect)]
      }, prevPos, curPos)
      
      
      ## Generate table of phased SNPs
      
      
      phasedList[[i]] <- table( unlist(matchedSnp.list) )
      i <- i + 1
    }
    
    prevPos <- curPos
  }
  return(phasedList)
}


c4.record_phasing <- function( hetPos.phased.list, sampleID, phasePath){
  cat(sampleID, file=phasePath, append=T)
  for( hetPos.table in hetPos.phased.list ){
    hetPos.vect <- names( hetPos.table )
    for(hetPos in hetPos.vect){
      hetDP <- as.integer( hetPos.table[hetPos] )
      cat(paste0('\t',hetPos,'_',hetDP), file=phasePath, append=T)
    }
  }
  cat('\n',file=phasePath,append=T)
}


c4.ins_pos <- function(startPos, cigarStr){
  
  output.insPos.list <- list()
  
  cigarOperationVect <- c('M','I','D')
  
  ## Format the start position
  startPos <- as.numeric(startPos)
  
  ## Split the cigar string up into a vector
  cigarVect <- strsplit(cigarStr,'')[[1]]
  
  ## Determine which position of the vector denote cigar operators
  operationPosVect <- which(cigarVect %in% cigarOperationVect)
  
  ## Initialize a list for storing the different operations and their positions
  operationList <- list()
  prevPos <- 1
  countInt <- 1
  
  ## Convert the cigar vect into the operation list
  for(operationPos in operationPosVect){
    operationList[[countInt]] <- cigarVect[prevPos:operationPos]
    prevPos <- operationPos+1
    countInt <- countInt+1
  }
  
  ## Initialize a start position for CIGAR operations
  output.iter <- 1
  prevPos <- 1
  subStrVect <- c()
  for(operationVect in operationList){
    
    ## Pull out the end position of the current operation
    operationEndPos <- as.integer(paste0(operationVect[1:(length(operationVect)-1)], collapse='')) + prevPos - 1
    
    ## Pull out the operation type
    operationTypeChar <- operationVect[length(operationVect)]
    
    ## If the end position is not an integer something went horribly wrong
    if(is.na(operationEndPos) | !(operationTypeChar %in% cigarOperationVect)){
      stop('c4.ins_pos something weird happened for ',readName)
    }
    
    if(operationTypeChar == 'M'){
      prevPos <- operationEndPos + 1
    }else if(operationTypeChar == 'I'){
      output.insPos.list[[output.iter]] <- prevPos:operationEndPos + startPos
      output.iter <- output.iter + 1
      
      ## Otherwise we ignore it
      prevPos <- operationEndPos + 1
    }else if(operationTypeChar == 'D'){
      
    }
  }
  
  return(output.insPos.list)
  
}
c4.exon_ins_pos_list <- function(samTable, c4.exonCoordVect){
  
  first.insIndex <- grepl('I', samTable$cigar_string)
  
  exon.insPos.list <- list()
  temp.iter <- 1
  
  first.insPos.list <- lapply( which( first.insIndex ) , function(x){
    x.cigarStr <- samTable$cigar_string[x]
    x.startPos <- samTable$startPos[x]
    
    x.output <- c4.ins_pos( x.startPos, x.cigarStr )
    
    x.isExon <- which( unlist( lapply( x.output, function(y){
      any( y %in% c4.exonCoordVect )
    }) ) )
    
    x.isExon.bool <- length(x.isExon) > 0
    
    if(x.isExon.bool){
      exon.insPos.list[[temp.iter]] <<- unlist( lapply( x.output[x.isExon], function(y){
        paste0(y, collapse='.')
      }) )
      
      temp.iter <<- temp.iter + 1
    }
    
    return(x.output)
  })
  
  return(exon.insPos.list)
}
c4.record_exon_ins <- function( total.exonINS.list, resultsDirectory, sampleID, indelPath){
  
  if( length(total.exonINS.list) > 0 ){
    exonINS.table <- table(unlist(total.exonINS.list))
    entryVect <- names( exonINS.table )
    
    for(insPos in entryVect){
      posDP <- exonINS.table[[insPos]]
      textStr <- paste0(sampleID,'\t',insPos,'\t',posDP,'\tINS\n')
      cat(textStr, file=indelPath, append=TRUE)
    }
  }
}
