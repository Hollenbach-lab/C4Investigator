
### Check to make sure bowtie2is accessible
bowtie2 <- system2('which', c('bowtie2'), stdout=T, stderr=T)
check.system2_output(bowtie2, 'bowtie2 not found')

## Check to make sure samtools is accessible
samtools <- system2('which', c('samtools'), stdout=T, stderr=T)
check.system2_output(samtools, 'samtools not found')

## This function runs a bowtie2 paired-end alignment
bowtie2.default <- function(bowtie2_command, reference_index, threads, currentSample, resultsDirectory){
  ## Intitialize an output path for the SAM file
  currentSample$samPath <- file.path(resultsDirectory,paste0(currentSample$name,'.sam'))
  
  ## Building up the run command
  optionsCommand <- c(paste0('-x ',reference_index),
                      '-I 0', '-X 1000', paste0('-p ',threads),
                      paste0('-1 ',currentSample$fastq1path),
                      paste0('-2 ',currentSample$fastq2path),
                      paste0('-S ',currentSample$samPath))
  
  ## Run the bowtie2 alignment command
  cat('\n\n',bowtie2_command,optionsCommand)
  output.sampleAlign <- system2(bowtie2_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.sampleAlign, 'bowtie2 gc alignment failed')
  
  ## Print the bowtie2 output
  cat('\n',paste0(output.sampleAlign, collapse='\n'))
  
  ## Check to make sure the SAM file actually exists
  currentSample$samPath <- normalizePath(currentSample$samPath, mustWork=T)
  
  cat('\n\nSuccessfully aligned',currentSample$name,'to',reference_index)
  
  return(output.sampleAlign)
}

## This function checks sam files for headers
check.samfile_header_exists <- function(currentSample){
  if(!file.exists(currentSample$samPath)){
    stop('This sam file does not exist')
  }
  
  con = file(currentSample$samPath, "r")
  
  line = readLines(con, n = 1)
  
  close(con)
  
  return(strsplit(line,'\t')[[1]][1] == '@SQ')
}

## This function counts the number of header lines in a SAM file (using '@SQ')
samfile.count_header_lines <- function(currentSample){
  if(!file.exists(currentSample$samPath)){
    stop('This sam file does not exist')
  }
  
  headerLines <- 0
  con = file(currentSample$samPath, "r")
  
  while(TRUE){
    
    line = readLines(con, n = 1)
    #if(strsplit(line,'\t')[[1]][1] != "@SQ"){
    #  break
    #}
    
    if(substr(line,1,1) != '@'){
      break
    }
    
    headerLines <- headerLines + 1
  }
  close(con)
  
  return(headerLines)
}

## This function generates a header for headerless sam files
samtools.generate_header <- function(samtools_command, currentSample, resultsDirectory, threads, referenceIndex){
  
  ## Initialize the output path for the temp sam file
  currentSample[['tempSam']] <- file.path(resultsDirectory,paste0(currentSample$name,'.header.sam'))
  
  ## Build up the run command
  optionsCommand <- c('view',paste0('-@',threads),'-ht',referenceIndex,
                      currentSample$samPath, '-o', currentSample$tempSam)
  
  cat('\n\n',samtools_command, optionsCommand)
  output.samGenHead <- system2(samtools_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.samGenHead, 'samtools header generation failed')
  
  ## Print the output
  cat('\n',paste0(output.samGenHead),collapse='\n')
  
  ## Make sure the temp file actually exists
  currentSample[['tempSam']] <- normalizePath(currentSample$tempSam,mustWork=T)
  
  ## Remove the headerless sam file
  file.remove(currentSample$samPath)
  
  ## Rename the header sam file to the normal samPath
  file.rename(currentSample$tempSam, currentSample$samPath)
  
  cat('\n\nSuccessfully generated a header for',currentSample$samPath)
  
  return(currentSample)
}

## This function runs a samtools sam to bam conversion
samtools.sam_to_bam <- function(samtools_command, currentSample, resultsDirectory, threads){
  
  ## Initialize an output path for the BAM file
  currentSample[['bamPath']] <- file.path(resultsDirectory,paste0(currentSample$name,'.bam'))
  
  ## Building up the run command
  optionsCommand <- c('view',paste0('-@', threads),
                      currentSample$samPath, '-o', currentSample$bamPath)
  
  cat('\n\n',samtools_command, optionsCommand)
  output.bamConv <- system2(samtools_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.bamConv, 'samtools sam to bam conversion failed')
  
  ## Print the conversion output
  cat('\n',paste0(output.bamConv), collapse='\n')
  
  ## Check to make sure the BAM file actually exists
  currentSample[['bamPath']] <- normalizePath(currentSample$bamPath, mustWork=T)
  
  cat('\n\nSuccessfully converted',currentSample$samPath,'to',currentSample$bamPath)
  
  return(currentSample)
}


## This function runs a samtools bam to sam conversion
samtools.bam_to_sam <- function(samtools_command, currentSample, resultsDirectory, threads){
  
  ## Building up the run command
  optionsCommand <- c('view',paste0('-@', threads),'-h',
                      currentSample$bamPath, '-o', currentSample$samPath)
  
  cat('\n\n',samtools_command, optionsCommand)
  output.samConv <- system2(samtools_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.samConv, 'samtools bam to sam conversion failed')
  
  ## Print the conversion output
  cat('\n',paste0(output.samConv), collapse='\n')
  
  ## Check to make sure the BAM file actually exists
  currentSample[['samPath']] <- normalizePath(currentSample$samPath, mustWork=T)
  
  cat('\n\nSuccessfully converted',currentSample$bamPath,'to',currentSample$samPath)
  
  return(currentSample)
}

## This function sorts a bam file
samtools.sort <- function(samtools_command, currentSample, resultsDirectory, threads){
  ## Initialize an output path for the BAM file
  currentSample[['sortedBamPath']] <- file.path(resultsDirectory,paste0(currentSample$name,'.sorted.bam'))
  
  ## Building up the run command
  optionsCommand <- c('sort',paste0('-@', threads),
                      currentSample$bamPath, '-o', currentSample$sortedBamPath)
  
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

## This function sorts a bam file
samtools.name_sort <- function(samtools_command, currentSample, resultsDirectory, threads){
  ## Initialize an output path for the BAM file
  currentSample[['nameSortedBamPath']] <- file.path(resultsDirectory,paste0(currentSample$name,'.nameSorted.bam'))
  
  ## Building up the run command
  optionsCommand <- c('sort',paste0('-@', threads),'-n',
                      currentSample$bamPath, '-o', currentSample$nameSortedBamPath)
  
  cat('\n\n',samtools_command, optionsCommand)
  output.bamSort <- system2(samtools_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.bamSort, 'samtools BAM sorting failed')
  
  ## Print the conversion output
  cat('\n',paste0(output.bamSort), collapse='\n')
  
  ## Check to make sure the sorted BAM file exists
  currentSample[['nameSortedBamPath']] <- normalizePath(currentSample$nameSortedBamPath, mustWork=T)
  
  cat('\n\nSuccessfully sorted',currentSample$nameSortedBamPath)
  
  return(currentSample)
}

# This function sorts a fixmate bam file
samtools.coord_sort <- function(samtools_command, currentSample, resultsDirectory, threads){
  ## Initialize an output path for the BAM file
  currentSample[['sortedBamPath']] <- file.path(resultsDirectory,paste0(currentSample$name,'.sorted.bam'))
  
  ## Building up the run command
  optionsCommand <- c('sort',paste0('-@', threads),
                      currentSample$fixmatePath, '-o', currentSample$sortedBamPath)
  
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

## This function fixes mates
samtools.fixmate <- function(samtools_command, currentSample, resultsDirectory){
  ## Initialize an output path for the BAM file
  currentSample[['fixmatePath']] <- file.path(resultsDirectory,paste0(currentSample$name,'.fixmate.bam'))
  
  ## Building up the run command
  optionsCommand <- c('fixmate','-m',
                      currentSample$nameSortedBamPath, currentSample$fixmatePath)
  
  cat('\n\n',samtools_command, optionsCommand)
  output.fixmate <- system2(samtools_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.fixmate, 'samtools fixmate failed')
  
  ## Print the conversion output
  cat('\n',paste0(output.fixmate), collapse='\n')
  
  ## Check to make sure the sorted BAM file exists
  currentSample[['fixmatePath']] <- normalizePath(currentSample$fixmatePath, mustWork=T)
  
  cat('\n\nSuccessfully sorted',currentSample$fixmatePath)
  
  return(currentSample)
}

## This function marks duplicate reads
samtools.markdup <- function(samtools_command, currentSample, resultsDirectory){
  ## Initialize an output path for the BAM file
  currentSample[['markdupPath']] <- file.path(resultsDirectory,paste0(currentSample$name,'.markdup.bam'))
  
  ## Building up the run command
  optionsCommand <- c('markdup','-s',
                      currentSample$sortedBamPath, currentSample$markdupPath)
  
  cat('\n\n',samtools_command, optionsCommand)
  output.markdup <- system2(samtools_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.markdup, 'samtools fixmate failed')
  
  ## Print the conversion output
  cat('\n',paste0(output.markdup), collapse='\n')
  
  ## Check to make sure the sorted BAM file exists
  currentSample[['markdupPath']] <- normalizePath(currentSample$markdupPath, mustWork=T)
  
  cat('\n\nSuccessfully markdup',currentSample$markdupPath)
  
  return(currentSample)
}

## This function indexes a sorted BAM file
samtools.index <- function(samtools_command, currentSample, resultsDirectory){
  ## Initialize an output path for the BAM file
  currentSample[['indexBamPath']] <- file.path(resultsDirectory,paste0(currentSample$name,'.sorted.bai'))
  
  ## Building up the run command
  optionsCommand <- c('index',currentSample$sortedBamPath,currentSample$indexBamPath)
  
  cat('\n\n',samtools_command, optionsCommand)
  output.bamIndex <- system2(samtools_command, optionsCommand, stdout=T, stderr=T)
  check.system2_output(output.bamIndex, 'samtools BAM indexing failed')
  
  ## Print the conversion output
  cat('\n',paste0(output.bamIndex), collapse='\n')
  
  ## Check to make sure the sorted BAM file exists
  currentSample[['indexBamPath']] <- normalizePath(currentSample$indexBamPath, mustWork=T)
  
  cat('\n\nSuccessfully indexed',currentSample$indexBamPath)
  
  return(currentSample)
}

## This function calulates the sequence depth along defined coordinates
samtools.depth <- function(samtools_command, currentSample, resultsDirectory, regionCoordinates){
  ## Initialize an output path for the BAM file
  currentSample[['depthPath']] <- file.path(resultsDirectory,paste0(currentSample$name,'.depth'))
  
  ## Building up the run command
  optionsCommand <- c('depth',paste0('-r ',regionCoordinates),'-aa',currentSample$sortedBamPath)
  
  cat('\n\n',samtools_command, optionsCommand)
  output.bamDepth <- system2(samtools_command, optionsCommand, stdout=currentSample$depthPath, stderr=currentSample$depthPath)
  check.system2_output(output.bamDepth, 'samtools BAM indexing failed')
  
  ## Check to make sure the sorted BAM file exists
  currentSample[['depthPath']] <- normalizePath(currentSample$depthPath, mustWork=T)
  
  cat('\n\nSuccessfully calculated depth of',currentSample$depthPath)
  
  return(currentSample)
}

## This function calulates the sequence depth for all reference positions
samtools.all_depth <- function(samtools_command, currentSample, resultsDirectory){
  ## Initialize an output path for the BAM file
  currentSample[['allDepthPath']] <- file.path(resultsDirectory,paste0(currentSample$name,'.allDepth'))
  
  ## Building up the run command
  optionsCommand <- c('depth','-a',currentSample$sortedBamPath)
  
  cat('\n\n',samtools_command, optionsCommand)
  output.bamDepth <- system2(samtools_command, optionsCommand, stdout=currentSample$allDepthPath, stderr=currentSample$allDepthPath)
  check.system2_output(output.bamDepth, 'samtools BAM indexing failed')
  
  ## Check to make sure the sorted BAM file exists
  currentSample[['allDepthPath']] <- normalizePath(currentSample$allDepthPath, mustWork=T)
  
  cat('\n\nSuccessfully calculated depth of',currentSample$allDepthPath)
  
  return(currentSample)
}

## This function calculates sequence depth along bedfile coordinates
samtools.bed_depth <- function(samtools_command, currentSample, resultsDirectory, bedPath){
  ## Initialize an output path for the BAM file
  currentSample[['depthPath']] <- file.path(resultsDirectory,paste0(currentSample$name,'.depth'))
  
  ## Building up the run command
  optionsCommand <- c('depth',paste0('-b ',bedPath),'-aa',currentSample$sortedBamPath)
  
  cat('\n\n',samtools_command, optionsCommand)
  output.bamDepth <- system2(samtools_command, optionsCommand, stdout=currentSample$depthPath, stderr=currentSample$depthPath)
  check.system2_output(output.bamDepth, 'samtools BAM indexing failed')
  
  ## Check to make sure the sorted BAM file exists
  currentSample[['depthPath']] <- normalizePath(currentSample$depthPath, mustWork=T)
  
  cat('\n\nSuccessfully calculated depth of',currentSample$depthPath)
  
  return(currentSample)
}

## This function reads in a SAM file with header to a data.table (header rows must be skipped)
samfile.read_whole_genome_sam <- function(sam_path, rows_to_skip=26){
  cat('\nReading in the SAM file.')
  
  ## Make sure the SAM file can be read in
  sam_path <- normalizePath(sam_path, mustWork=T)
  
  ### rows_to_skip should include all header rows
  #output.samTable <- fread(sam_path, sep='\t', stringsAsFactors=F, check.names=F, fill=T, skip=rows_to_skip)
  
  #output.samTable <- read.table(sam_path, sep='\t', col.names=1:25, stringsAsFactors=F, check.names=F, fill=T, skip=rows_to_skip, nrows = 30000000)
  output.samTable <- read.table(sam_path, sep='\t', stringsAsFactors=F, check.names=F, fill=T, skip=rows_to_skip,col.names=1:25, quote="")
  
  ## Convert the dataframe to a datatable for faster access
  output.samTable <- as.data.table(output.samTable)
  
  ## Name the columns that are used for downstream analysis
  colnames(output.samTable)[1] <- 'read_name'
  colnames(output.samTable)[2] <- 'sam_flag'
  colnames(output.samTable)[3] <- 'reference_name'
  colnames(output.samTable)[4] <- 'ref_pos'
  colnames(output.samTable)[6] <- 'cigar_string'
  colnames(output.samTable)[10] <- 'read_seq'
  
  ## Convert the alignment scores into integers
  alignmentScoreList <- as.integer(tstrsplit(output.samTable$'12', ':', fixed=TRUE)[[3]])
  
  ## Save alignment scores to their own columns
  cat('\n\tSetting alignment scores.')
  output.samTable[,alignment_score := alignmentScoreList]
  
  ## Use the reference name conversion list to apply common gene names
  output.samTable$locus <- unlist(lapply(output.samTable$reference_name, function(x){ referenceKeyList[[ as.character(x) ]] }))
  
  ## Pull out an index of reads that aligned to C4
  C4ReadsPosVect <- grep('C4',output.samTable$locus,fixed=T)
  
  ## Initialize a column of the data table for storing whether the read aligned to C4 or not
  output.samTable$isC4 <- F
  output.samTable$isC4[C4ReadsPosVect] <- T ## Set the C4 bool to T for all C4 indexed reads
  
  return(output.samTable)
}

## This function interprets SAM hex flags
samtable.flag_convert <- function(samFlag){
  flagList <- list(readPaired=F,properPairMapping=F,readUnmapped=F,mateUnmapped=F,
                   readReverseStrand=F,mateReverseStrand=F,firstInPair=F,secondInPair=F,
                   notPrimaryAlignment=F,readFailsQualityChecks=F,readIsPcrOrOpticalDuplicate=F,
                   supplementaryAlignment=F)
  
  #read paired (0x1)
  #read mapped in proper pair (0x2)
  #read unmapped (0x4)
  #mate unmapped (0x8)
  #read reverse strand (0x10)
  #mate reverse strand (0x20)
  #first in pair (0x40)
  #second in pair (0x80)
  #not primary alignment (0x100)
  #read fails platform/vendor quality checks (0x200)
  #read is PCR or optical duplicate (0x400)
  #supplementary alignment (0x800)
  
  
  hexInt <- as.integer(samFlag)
  
  if(hexInt >= 2048){
    hexInt <- hexInt - 2048
    flagList$supplementaryAlignment <- T
  }
  
  if(hexInt >= 1024){
    hexInt <- hexInt - 1024
    flagList$readIsPcrOrOpticalDuplicate <- T
  }
  
  if(hexInt >= 512){
    hexInt <- hexInt - 512
    flagList$readFailsQualityChecks <- T
  }
  
  if(hexInt >= 256){
    hexInt <- hexInt - 256
    flagList$notPrimaryAlignment <- T
  }
  
  if(hexInt >= 128){
    hexInt <- hexInt - 128
    flagList$secondInPair <- T
  }
  
  if(hexInt >= 64){
    hexInt <- hexInt - 64
    flagList$firstInPair <- T
  }
  
  if(hexInt >= 32){
    hexInt <- hexInt - 32
    flagList$mateReverseStrand <- T
  }
  
  if(hexInt >= 16){
    hexInt <- hexInt - 16
    flagList$readReverseStrand <- T
  }
  
  if(hexInt >= 8){
    hexInt <- hexInt - 8
    flagList$mateUnmapped <- T
  }
  
  if(hexInt >= 4){
    hexInt <- hexInt - 4
    flagList$readUnmapped <- T
  }
  
  if(hexInt >= 2){
    hexInt <- hexInt - 2
    flagList$properPairMapping <- T
  }
  
  if(hexInt >= 1){
    hexInt <- hexInt - 1
    flagList$readPaired <- T
  }
  
  return(flagList)
}
