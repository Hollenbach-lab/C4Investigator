
kirLocusList <- c('KIR2DL1','KIR2DL2','KIR2DL3','KIR2DL4','KIR2DL5A','KIR2DL5B',
                  'KIR2DS1','KIR2DS2','KIR2DS3','KIR2DS4','KIR2DS5','KIR2DP1',
                  'KIR3DL1','KIR3DL2','KIR3DL3','KIR3DS1','KIR3DP1')

nucListConv <- list('A'=1,
                    'T'=2,
                    'C'=3,
                    'G'=4,
                    '.'=5)

## Creating the results direcotry object class
pathObj <- setRefClass("pathObj",
                       fields=list(name='character',
                                   path='character',
                                   subDir='list'))

pathObj$methods(
  dirGen= function(){
    if(!file.exists(path)) dir.create(path, recursive = T, showWarnings = F)
    path <<- normalizePath(path, mustWork = T)
    invisible(path)
  }
)

## This function checks to make sure the output of system2 is valid
check.system2_output <- function(system2_output, system2_error){
  
  ## Checking if the attributes of system2_command are NULL, if not the command was not found
  if(!is.null(attributes(system2_output))){
    cat('\n',system2_output)
    stop(system2_error, '. Stopping program.')
  }
}

general.identify_uppercase <- function(str){
  return( str == toupper(str) )
}

## This function finds samples in sampleDirectory and turns them into objects for downstream use
sequence.paired_sample_objects <- function(sample_directory, fastq_pattern='fastq', resultsDirectory){
  #####
  ## This function takes in a directory and a file name pattern and attempts to pair fastq files
  ## Returns a list of sample objects that contain the paired fastq file names
  #####
  
  cat("\nAttempting automatic fastq pairing in", sample_directory, "using", fastq_pattern)
  
  ## Find all the files in sampleDirectory that match fastqPattern
  unpairedFastqList <- list.files(path=sample_directory, pattern=fastq_pattern)
  
  ## To pair reads, we will split the file names by fastqPattern, then continuously chop a
  ## character off the end of each name until the number of unique names is exactly half of the total names
  
  ## Setting up an initial fastq list that splits the files names by fastqPattern
  strList <- sapply(unpairedFastqList, function(x) str_split(x, fastq_pattern)[[1]][1])
  
  ## Setting the maximum number of times to chop off the last character to the length of the shortest name
  maxChop <- min(sapply((strList), nchar))
  
  ## Iterate from 0 to maxChop, cutting off i characters from the file names each time
  for(i in 0:maxChop){
    
    ## In each iteration, the file names are reset. There is no particular reason I implemented it this way
    subStrList <- strList
    
    ## Cut off i characters from the end of each fastq file name
    subStrList <- sapply(subStrList, function(x) substr(x, 1, nchar(x)-i))
    
    ## After cutting, determine the unique names
    uniqueFastqList <- unique(subStrList)
    
    ## If the number of unique names is exactly half of the total names, then cutting should be finished
    ## and it is time to move on to matching!
    if(length(uniqueFastqList) == (length(subStrList)/2)){
      break
    }
    
    ## Pairing failed if i reaches maxChop. This will raise a fatal error.
    if(i == maxChop){
      stop("Was not able to pair fastq file names, please check that fastqPattern and sampleDirectory are set correctly.")
    }
  }
  
  ## Initialize a list for storing the paired fastq file names
  pairedFastqList <- list()
  
  ## Iterate through the unique fastq names to make pairings
  for(fastqName in uniqueFastqList){
    
    ## Pull out the matches for fastqName in the subStrList
    fastqMatches <- subStrList[fastqName == subStrList]
    
    ## Determine how many matches there are for fastqName in the subStrList
    matchCount <- length(fastqMatches)
    
    ## Stop the program if more or less than 2 matches are found for fastqName
    if(matchCount != 2){
      cat('\n',names(fastqMatches))
      stop('Auto fastq matching failed due to an improper number of matches for ',fastqName)
    }
    
    ## Save the file names in the pairedFastqList under the unique name
    pairedFastqList[[fastqName]] <- names(fastqMatches)
  }
  
  cat("\nFound", length(uniqueFastqList), "samples in", sample_directory)
  
  ## Creating the sample object class
  sample <- setRefClass("sample",
                        fields=list(name='character',
                                    fastq1path='character',
                                    fastq2path='character',
                                    gzip='logical',
                                    samPath='character',
                                    bamPath='character'))
  
  ## Initializing a sample object list. This will be returned
  output.sampleList <- list()
  
  for(i in 1:length(pairedFastqList)){
    
    ## Pulling the current working element out of the list
    pairedFastq <- pairedFastqList[i]
    
    ## Creating an absolute path to the first fastq file
    fastq1path <- normalizePath(file.path(sample_directory, pairedFastq[[1]][1]), mustWork=T)
    
    ## Creating a absolute path to the second fastq file
    fastq2path <- normalizePath(file.path(sample_directory, pairedFastq[[1]][2]), mustWork=T)
    
    ## Checking if the first fastq file is gzipped
    gzip <- substr(fastq1path, nchar(fastq1path)-2, nchar(fastq1path)) == '.gz'
    
    ## Fill in the path to the alignment file (it may or may not be present)
    samPath <- file.path(resultsDirectory,paste0(names(pairedFastq),'.sam'))
    
    ## Fill in the path to the alignment file (it may or may not be present)
    bamPath <- file.path(resultsDirectory,paste0(names(pairedFastq),'.bam'))
    
    ## Building a sample object and adding it to sampleList
    output.sampleList[[names(pairedFastq)]] <- sample(name=names(pairedFastq),fastq1path=fastq1path,fastq2path=fastq2path,gzip=gzip,samPath=samPath,bamPath=bamPath)
  }
  
  cat("\nAll samples were successfully paired")
  return(output.sampleList)
}



