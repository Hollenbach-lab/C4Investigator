# Copy computation --------------------------------------------------------
c4.calc_total_copy <- function(median_mhc_c4, median_tnxb){
  dum <- (median_mhc_c4 / median_tnxb)*4
  c4_copy <- dum - (dum*0.1)
  return(c4_copy)
}

c4.calc_AB_copy <- function(c4_copy, c4a_g1, c4a_g2, c4b_g1, c4b_g2){
  
  rouded_copy <- round(c4_copy, digits=0)
  
  c4a_g1_copy <- (c4a_g1 / (c4a_g1+c4b_g1))*rouded_copy
  c4b_g1_copy <- (c4b_g1 / (c4a_g1+c4b_g1))*rouded_copy
  
  c4a_g2_copy <- (c4a_g2 / (c4a_g2+c4b_g2))*rouded_copy
  c4b_g2_copy <- (c4b_g2 / (c4a_g2+c4b_g2))*rouded_copy
  
  #c4a_copy <- (c4a_g1_copy + c4a_g2_copy)/2
  #c4b_copy <- (c4b_g1_copy + c4b_g2_copy)/2
  c4a_copy <- c4a_g1_copy
  c4b_copy <- c4b_g1_copy
  
  cat('\n\nC4A:',
      round(c4a_copy,2),
      paste0('\n\tg1:',round(c4a_g1_copy,2))
      #paste0('\n\tg2:',round(c4a_g2_copy,2))
  )
  
  cat('\n\nC4B:',
      round(c4b_copy,2),
      paste0('\n\tg1:',round(c4b_g1_copy,2))
      #paste0('\n\tg2:',round(c4b_g2_copy,2))
  )
  
  
  return(list('c4a_g1_copy'=round(c4a_g1_copy,0), 
              'c4rg_g2_copy'=round(c4a_g2_copy,0), 
              'c4b_g1_copy'=round(c4b_g1_copy,0), 
              'c4ch_g2_copy'=round(c4b_g2_copy,0)))
  
}

c4.calc_SL_copy <- function(c4_copy, med_ins_depth, med_del_depth){
  rounded_copy <- round(c4_copy, digits=0)
  
  c4s_copy <- round( (med_del_depth/(med_ins_depth+med_del_depth))*rounded_copy, digits=0)
  c4l_copy <- round( (med_ins_depth/(med_ins_depth+med_del_depth))*rounded_copy, digits=0)
  
  #c4s_copy <- min(round((med_del_depth/med_c4_depth)*rounded_copy,digits=0), rounded_copy)
  #c4l_copy <- rounded_copy - c4s_copy
  return(list('c4s'=c4s_copy,
              'c4l'=c4l_copy))
}

c4.separate_recombinants <- function( c4DT ){
  
  recombBoolVect <- c( round(c4DT$C4Ag1_copy, 0) != round(c4DT$C4Ag2_copy, 0)) | (round(c4DT$C4Bg1_copy, 0) != round(c4DT$C4Bg2_copy, 0) )
  
  recombDT <- c4DT[recombBoolVect,]
  otherDT <- c4DT[!recombBoolVect,]
  
  return(list('rcDT'=recombDT,
              'normDT'=otherDT))
}

## Returns boolean
is_nuc <- function(chr){
  nuc_list <- c('A', 'T', 'C', 'G', '.')
  return(as.character(chr) %in% nuc_list)
}

## Counts how many unique nucleotides are in a vector
num_unique_nuc <- function(vector_to_check){
  return(sum(is_nuc(names(table(vector_to_check)))))
}

c4.generate_feature_list <- function( ref.fa ){
  if( !file.exists(ref.fa) ){
    stop(paste('reference file does not exist here:',ref.fa))
  }
  ref.dt <- fread(file=ref.fa,header=F)
  ref.seq <- ref.dt[2,][[1]]
  
  if( nchar(ref.seq) != 21058 ){
    stop(paste('reference sequence is not expected length'))
  }
  
  ref.seq.vect <- strsplit(ref.seq,'')[[1]]
  exonCoord.index <- which( ref.seq.vect == toupper(ref.seq.vect) )
  
  boundary.vect <- which( diff(exonCoord.index) > 1 )
  
  initial.pos <- 285
  exon.number <- 1
  intron.number <- 0
  geneFeature.list <- list()
  for( boundary in boundary.vect ){
    end.pos <- exonCoord.index[boundary]
    
    ref.seq.vect[(initial.pos-1):(end.pos+1)]
    
    feat.index <- initial.pos:end.pos
    
    if(intron.number == 0){
      geneFeature.list[['5UTR']] <- paste0('5UTR_',1:(initial.pos-1))
      geneFeature.list[[paste0('E',exon.number)]] <- paste0('E',exon.number,'_',1:length( feat.index ))
    }else{
      geneFeature.list[[paste0('I',intron.number)]] <- paste0('I',intron.number,'_', 1:length( (length(unlist(geneFeature.list))+1):(initial.pos-1) ))
      geneFeature.list[[paste0('E',exon.number)]] <- paste0('E',exon.number,'_',1:(length( feat.index )))
    }
    
    if( boundary == boundary.vect[length(boundary.vect)] ){
      geneFeature.list[['I40']] <- paste0('I40_',1:length((feat.index[length(feat.index)]+1):(exonCoord.index[(boundary+1)]-1)))
      geneFeature.list[['E41']] <- paste0('E41_',1:length( exonCoord.index[(boundary+1)]:20717 ))
      geneFeature.list[['3UTR']] <- paste0('3UTR_',1:length(20718:(length(ref.seq.vect))))
    }
    
    exon.number <- exon.number + 1
    intron.number <- intron.number + 1
    initial.pos <- exonCoord.index[boundary+1]
  }
  
  geneFeature.vect <- unlist(geneFeature.list)
  
  #ref.seq.vect[ grep('5UTR_',geneFeature.vect,fixed=T) ]
  #ref.seq.vect[ grep('E1_',geneFeature.vect,fixed=T) ]
  #ref.seq.vect[ grep('I1_',geneFeature.vect,fixed=T) ]
  #ref.seq.vect[ grep('E2_',geneFeature.vect,fixed=T) ]
  
  exonErrors <- which( ref.seq.vect[ grep('E',geneFeature.vect) ] != toupper(ref.seq.vect[ grep('E',geneFeature.vect) ]) )
  intronErrors <- which( ref.seq.vect[ grep('I',geneFeature.vect) ] != tolower(ref.seq.vect[ grep('I',geneFeature.vect) ]) )
  
  if(length(exonErrors) > 0 ){
    stop('non exonic positions in exon index')
  }
  
  if(length(intronErrors) > 0 ){
    stop('non intronic positions in intron index')
  }
  
  if(length(unlist(geneFeature.list)) != 21058 ){
    stop('gene feature list is not expected length')
  }
  
  return( geneFeature.list )
}
c4.feature.list <- c4.generate_feature_list('resources/c4only_onelines_oneDel_bShort.fasta')
c4.feature.vect <- 1:length(unlist(c4.feature.list))
names(c4.feature.vect) <- unlist(c4.feature.list)
ins.feature.vect <- c4.feature.vect[2860:9230]

important.hetPos.vect <- c("E25_64","E26_55","E26_129","E26_132","E26_140","E26_143","E26_145","E28_23","E28_111","E28_116","E28_125","E28_126","I28_14","I28_19")

nucListConv <- list('A'=1,
                    'T'=2,
                    'C'=3,
                    'G'=4,
                    '.'=5)

c4.read_dp_csv <- function(sampleID, rawData.dir, c4.feature.list, nucListConv, onlyExon=T){
  sample.dir <- list.files(path=file.path(rawData.dir,'home'), pattern=sampleID,recursive = T,include.dirs = T,full.names = T)[[1]]
  sample.dpDF.path <- list.files(path=sample.dir,pattern = 'c4_dp.csv',recursive=T,full.names=T)
  
  sample.dt <- fread(file=sample.dpDF.path,header=T,nrows=5)
  sample.dt <- sample.dt[,2:ncol(sample.dt)]
  
  colnames(sample.dt) <- unlist(c4.feature.list)
  if(onlyExon){
    exon.index <- as.vector( grep('E',unlist(c4.feature.list),fixed=T,value=T) )
    sample.dt <- sample.dt[,exon.index,with=F]
  }
  #sample.dt$nuc <- names((nucListConv))
  return(sample.dt)
}

c4.set_snp_dt <- function( sample.dt, hetRatio, sampleCopy, cutPos.vect ){
  totalDP.vect <- apply( sample.dt, 2, sum)
  copyDP.vect <- totalDP.vect / sampleCopy
  
  snpBool.mat <- sample.dt > copyDP.vect*hetRatio
  
  if(length(cutPos.vect) > 0){
    cutPos.nameVect <- names(cutPos.vect)
    #snpBool.mat <- snpBool.mat[,!colnames(snpBool.mat) %in% cutPos.nameVect]
    snpBool.mat[,cutPos.nameVect] <- F
  }
  
  a.index <- which( snpBool.mat[1,] )
  t.index <- which( snpBool.mat[2,] )
  c.index <- which( snpBool.mat[3,] )
  g.index <- which( snpBool.mat[4,] )
  del.index <- which( snpBool.mat[5,] )
  
  currentSample.snp.dt <- as.data.table(matrix(data='',nrow=4,ncol=ncol(snpBool.mat)))
  colnames(currentSample.snp.dt) <- colnames(sample.dt)
  
  snp.tab <- table( c( names(a.index), names(t.index), names(c.index), names(g.index), names(del.index) ))
  nuc2.vect <- which( snp.tab == 2 )
  nuc3.vect <- which( snp.tab == 3 )
  nuc4.vect <- which( snp.tab == 4 )
  
  hetPos.vect <- c()
  snpNum.vect <- c(1)
  if( length(nuc2.vect) > 0 ){
    duoPos.vect <- names(nuc2.vect)
    hetPos.vect <- c(hetPos.vect, duoPos.vect)
    snpNum.vect <- c(snpNum.vect, 2)
  }
  
  if( length(nuc3.vect) > 0 ){
    triPos.vect <- names(nuc3.vect)
    hetPos.vect <- c(hetPos.vect, triPos.vect)
    snpNum.vect <- c(snpNum.vect, 3)
  }
  
  if( length(nuc4.vect) > 0 ){
    quadPos.vect <- names(nuc4.vect)
    hetPos.vect <- c(hetPos.vect, quadPos.vect)
    snpNum.vect <- c(snpNum.vect, 4)
  }
  
  snpNum <- 1
  cols <- setdiff( names(a.index), hetPos.vect )
  if( length(cols) > 0) set(currentSample.snp.dt, as.integer(snpNum), cols, 'A' )
  
  cols <- setdiff( names(t.index), hetPos.vect )
  if( length(cols) > 0) set(currentSample.snp.dt, as.integer(snpNum), cols, 'T' )
  
  cols <- setdiff( names(c.index), hetPos.vect )
  if( length(cols) > 0) set(currentSample.snp.dt, as.integer(snpNum), cols, 'C' )
  
  cols <- setdiff( names(g.index), hetPos.vect )
  if( length(cols) > 0) set(currentSample.snp.dt, as.integer(snpNum), cols, 'G' )
  
  cols <- setdiff( names(del.index), hetPos.vect )
  if( length(cols) > 0) set(currentSample.snp.dt, as.integer(snpNum), cols, '.' )
  
  if( length(nuc2.vect) > 0 ){
    for( pos in duoPos.vect ){
      pos.nucVect <- c()
      if( pos %in% names(a.index) ) pos.nucVect <- c(pos.nucVect, 'A')
      if( pos %in% names(t.index) ) pos.nucVect <- c(pos.nucVect, 'T')
      if( pos %in% names(c.index) ) pos.nucVect <- c(pos.nucVect, 'C')
      if( pos %in% names(g.index) ) pos.nucVect <- c(pos.nucVect, 'G')
      if( pos %in% names(del.index) ) pos.nucVect <- c(pos.nucVect, '.')
      set(currentSample.snp.dt, as.integer(1:2), pos, pos.nucVect)
    }
  }
  
  if( length(nuc3.vect) > 0 ){
    for( pos in triPos.vect ){
      pos.nucVect <- c()
      if( pos %in% names(a.index) ) pos.nucVect <- c(pos.nucVect, 'A')
      if( pos %in% names(t.index) ) pos.nucVect <- c(pos.nucVect, 'T')
      if( pos %in% names(c.index) ) pos.nucVect <- c(pos.nucVect, 'C')
      if( pos %in% names(g.index) ) pos.nucVect <- c(pos.nucVect, 'G')
      if( pos %in% names(del.index) ) pos.nucVect <- c(pos.nucVect, '.')
      set(currentSample.snp.dt, as.integer(1:3), pos, pos.nucVect)
    }
  }
  
  if( length(nuc4.vect) > 0 ){
    for( pos in quadPos.vect ){
      pos.nucVect <- c()
      if( pos %in% names(a.index) ) pos.nucVect <- c(pos.nucVect, 'A')
      if( pos %in% names(t.index) ) pos.nucVect <- c(pos.nucVect, 'T')
      if( pos %in% names(c.index) ) pos.nucVect <- c(pos.nucVect, 'C')
      if( pos %in% names(g.index) ) pos.nucVect <- c(pos.nucVect, 'G')
      if( pos %in% names(del.index) ) pos.nucVect <- c(pos.nucVect, '.')
      set(currentSample.snp.dt, as.integer(1:4), pos, pos.nucVect)
    }
  }
  
  return(currentSample.snp.dt)
}

c4.read_cov <- function(cov.path,covData.dir){
  dp.median <- read.table( file.path(covData.dir,cov.path), sep='\t' )[1,7]
  
  sampleID <- strsplit(cov.path,'_',fixed=T)[[1]][1]
  return(list('sampleID'=sampleID,'median'=dp.median))
}

processUniprot <- function(uniprotID){
  
  PTM.info <- GetProteinAnnontate(uniprotID,c('feature(CHAIN)', 'feature(CROSS LINK)', 'feature(DISULFIDE BOND)', 
                                              'feature(GLYCOSYLATION)', 'feature(INITIATOR METHIONINE)', 'feature(LIPIDATION)', 
                                              'feature(MODIFIED RESIDUE)', 'feature(PEPTIDE)', 'feature(PROPEPTIDE)',
                                              'feature(SIGNAL)', 'feature(TRANSIT)','feature(NATURAL VARIANT)'))
  out.list <- list()
  
  for( PTM.feature in names(PTM.info) ){
    
    if(PTM.feature == "Natural.variant"){
      featureStr <- gsub('; dbSNP:',' ;/dbSNP=', PTM.info[[PTM.feature]],fixed=T)
      featureStr <- gsub('(','',featureStr,fixed=T)
      featureStr <- gsub(')','',featureStr,fixed=T)
    }else if(PTM.feature == "Modified.residue"){
      featureStr <- gsub('; by', ' by', PTM.info[[PTM.feature]])
    }else{
      featureStr <- PTM.info[[PTM.feature]]
    }
    featureVect <- trimws(strsplit(featureStr,';',fixed=T)[[1]])
    
    if(length(featureVect) == 1 && featureVect == "NA") next
    
    feature.list <- list()
    for( feature.elem in featureVect ){
      
      if( !grepl('/',feature.elem,fixed=T) ){
        feature.id <- feature.elem
        feature.list[[feature.id]] <- list()
        
        feature.pos <- strsplit(feature.id,' ',fixed=T)[[1]][2]
        if( grepl('..',feature.pos,fixed=T) ){
          temp1 <- as.integer(strsplit(feature.pos,'..',fixed=T)[[1]])
          
          if( PTM.feature %in% c('Disulfide.bond','Cross.link')){
            feature.pos <- temp1
          }else{
            feature.pos <- temp1[1]:temp1[2]
          }
          
        }else{
          feature.pos <- as.integer(feature.pos)
        }
        
        feature.list[[feature.id]][['pos']] <- feature.pos
      }else{
        
        elem.input <- strsplit(strsplit(feature.elem,'/',fixed=T)[[1]][2], '=', fixed=T)[[1]]
        
        feature.list[[feature.id]][[elem.input[1]]] <- elem.input[2]
      }
    }
    
    out.list[[PTM.feature]] <- feature.list
  }
  
  PTM.info <- GetProteinAnnontate(uniprotID,c('sequence'))
  
  out.list[['Sequence']] <- PTM.info
  
  return(out.list)
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
      alleleName <- strsplit(currentLine, ' ',fixed=TRUE)[[1]][1]
      alleleName <- gsub('>','',alleleName,fixed = T)
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

general.identify_uppercase <- function(str){
  return( str == toupper(str) )
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

#### These sets of functions use the deletion index lists to convert read alignment coordinates to a common coordinate
## This function sets the universal read start position
run.setStartPos <- function(ref_name, ref_pos){
  return(deletionIndexList[[ref_name]][ref_pos])
}
## This function sets the universal read end position
run.setEndPos <- function(ref_name, ref_pos, currentReadLen){
  return(deletionIndexList[[ref_name]][ref_pos+currentReadLen-1])
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


## Format reads for snp phasing (fill in deletion positions with '.')
alleleSetup.read_formatter <- function( readName, currentSeq, refPos, currentLocus, cigarStr ){
  
  #cat('\n',readName,currentLocus,currentRef)
  cigarMod.del.bool <- grepl('D', cigarStr)
  cigarMod.ins.bool <- grepl('I', cigarStr)
  
  if( cigarMod.del.bool | cigarMod.ins.bool ){
    resultList <- samFormat.processCigarStr( cigarStr, currentSeq )
    currentSeq <- resultList$readSeq
    insIndex <- resultList$insIndex
  }
  
  currentDelIndex <- inverseDeletionIndexList[[currentLocus]]
  
  #cat('',length(currentDelIndex))
  
  seqLength <- nchar(currentSeq)
  delBool <- length(currentDelIndex) > 0
  
  delRegionIndexVect <- which( diff( currentDelIndex ) > 1 )
  
  startPos <- alleleSetup.del_start_offset( refPos, currentDelIndex, delRegionIndexVect )
  #cat('','passed startPos')
  ## Format the start position
  #startPos <- as.numeric(refPos)
  
  ## Format the end position
  endDelIndex <- currentDelIndex[ currentDelIndex > startPos ]
  delRegionIndexVect <- which( diff( endDelIndex ) > 1 )
  endPos <- startPos + alleleSetup.del_end_offset( startPos, seqLength, endDelIndex, delRegionIndexVect ) - 1
  
  ## Testing if there deletions found in the reference for the current read
  if( (endPos-startPos + 1) > seqLength ){
    
    delIndexList <- which(startPos:endPos %in% currentDelIndex)
    
    ## If the index list is empty, something went wrong
    if(length(delIndexList) == 0){
      stop('Something weird happened with the deletion index.')
    }
    
    ## Split apart the read sequence and add in deletion symbols
    for(delIndex in delIndexList){
      seqLength <- nchar(currentSeq)
      preString <- substr(currentSeq, 1, delIndex-1)
      postString <- substr(currentSeq, delIndex, seqLength)
      
      currentSeq <- paste0(preString, '.', postString)
    }
  }
  
  ## Create a full index that cooresponds to the sequence nucleotides and where they go
  fullIndex <- startPos:endPos
  
  ## Turn the sequence string into a list
  seqList <- strsplit(currentSeq,'')[[1]]
  
  if(length(seqList) != length(fullIndex)){
    cat('\n',startPos,endPos,locus,length(seqList),length(fullIndex))
  }
  
  names(seqList) <- fullIndex
  
  if( cigarMod.ins.bool ){
    for( insPos in names( insIndex ) ){
      insPosInt <- as.integer(insPos)
      #cat('\n',readName)
      # Mod to assign insertion sequence over subsequet del characters
      if( seqList[ (as.integer(insPos)+1) ] == '.' ){
        
        insSeqVect <- strsplit( insIndex[[insPos]], '')[[1]]
        insSeqVect <- insSeqVect[2:length(insSeqVect)]
        insSeqPosVect <- (insPosInt+1):(length(insSeqVect)+insPosInt)
        
        readDelIndex <- which( seqList == '.' )
        
        assign.index <- insSeqPosVect %in% readDelIndex
        
        if( any(assign.index) ){
          seqList[ insSeqPosVect[ assign.index ] ] <- insSeqVect[ assign.index ]
        }
        
        # If there is leftover insertion sequence, append it to the last assigned ins character
        if( any(!assign.index) ){
          lastInsPos <- insSeqPosVect[ assign.index ][ sum( assign.index ) ]
          seqList[ lastInsPos ] <- paste0( c(seqList[ lastInsPos ], insSeqVect[ !assign.index ]), collapse='')
        }
      }else{
        seqList[ insPosInt ] <- insIndex[[insPos]]
      }
    }
  }
  
  seqListList <- list( list(seqList) )
  return(seqListList)
}

cigarOperationVect <- c('M','I','D')
samFormat.processCigarStr <- function( cigarStr, readSeq ){
  # INDEL handling
  # deletions - add . character to read for deletion positions marked in CIGAR string
  # insertions - remove insertion characters from read, send read through the indexed vector formatting step, 
  #              then replace the insertion position with the insertion string ( will need to add the current
  #              read del index to correctly place )
  # reads with both - should work following the above steps
  
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
  prevPos <- 1
  outStr <- ''
  out.insIndex <- list()
  
  #operationVect <- operationList[[6]]
  for(operationVect in operationList){
    
    ## Pull out the end position of the current operation
    operationEndPos <- as.integer(paste0(operationVect[1:(length(operationVect)-1)], collapse='')) + prevPos - 1
    
    ## Pull out the operation type
    operationTypeChar <- operationVect[length(operationVect)]
    
    if(operationTypeChar == 'M'){
      prevPos <- operationEndPos + 1
    }else if(operationTypeChar == 'I'){
      
      preString <- substr(readSeq, 1, prevPos-1)
      insStr <- substr(readSeq, prevPos-1, operationEndPos)
      postString <- substr(readSeq, operationEndPos+1, nchar(readSeq) )
      
      out.insIndex[[as.character( prevPos-1 )]] <- insStr
      readSeq <- paste0( preString, postString )
      
      ## We make no changes to the position variables since the INS seq was removed
      #operationOffset <- operationOffset + nchar(insStr) - 1
      
      #prevPos <- operationEndPos + 1
    }else if(operationTypeChar == 'D'){
      
      seqLength <- nchar(readSeq)
      preString <- substr(readSeq, 1, prevPos-1)
      postString <- substr(readSeq, prevPos, seqLength)
      readSeq <- paste0(c( preString, rep('.', length( prevPos:operationEndPos )), postString), collapse='')
      
      prevPos <- operationEndPos + 1
      
    }else{
      stop('samFormat.adjustReadSeq something weird happened for ',readSeq)
    }
    
  }
  
  return(list('readSeq'=readSeq,'insIndex'=out.insIndex))
}


alleleSetup.del_start_offset <- function( refPos, currentDelIndex, delRegionIndexVect ){
  #cat('\n',refPos,'',length(currentDelIndex))
  
  if(length( delRegionIndexVect ) == 0 ){
    delRegionIndexVect <- length(currentDelIndex)
  }
  
  if(length(currentDelIndex) == 0){
    return(refPos)
  }else{
    
    #oldRefPos <- refPos
    delOffsetBool <- any( refPos >= currentDelIndex[1:delRegionIndexVect[1]] )
    #delOffset <- length(currentDelIndex)
    
    if( !delOffsetBool ){
      return(refPos)
    }
    
    delOffset <- length(currentDelIndex[1:delRegionIndexVect[1]])
    refPos <- refPos + delOffset
    
    if(length( currentDelIndex ) > delOffset ){
      currentDelIndex <- currentDelIndex[ (delOffset+1):length(currentDelIndex) ]
      delRegionIndexVect <- delRegionIndexVect - delOffset
      return( alleleSetup.del_start_offset( refPos, currentDelIndex, delRegionIndexVect[-1]) )
    }else{
      return(refPos)
    }
    
  }
  
}

alleleSetup.del_end_offset <- function( startPos, seqLength, endDelIndex, delRegionIndexVect ){
  #cat('\n',startPos,'',seqLength,'',length(endDelIndex))
  endPos <- startPos + seqLength
  
  if(length( delRegionIndexVect ) == 0 ){
    delRegionIndexVect <- length(endDelIndex)
  }
  
  if(length(endDelIndex) == 0){
    return( seqLength )
  }else{
    
    delOffsetBool <- any( endPos >= endDelIndex )
    
    if( !delOffsetBool ){
      return( seqLength )
    }
    
    delOffset <- length( endDelIndex[1:delRegionIndexVect[1]] )
    
    seqLength <- seqLength + delOffset
    
    if( length(endDelIndex) > delOffset ){
      endDelIndex <- endDelIndex[ (delOffset+1):length(endDelIndex) ]
      delRegionIndexVect <- delRegionIndexVect - delOffset
      return( alleleSetup.del_end_offset( startPos, seqLength, endDelIndex, delRegionIndexVect[-1] ) )
    }else{
      return( seqLength )
    }
  }
}

c4.format_full_readTable <- function( readTable.list ){
  
  readTable.posVect <- as.integer(names(unlist(readTable.list))) 
  full.readTable.posVect <- min(readTable.posVect):max(readTable.posVect)
  names(full.readTable.posVect) <- as.character(full.readTable.posVect)
  
  empty.posVect <- setdiff(names(full.readTable.posVect),as.character(readTable.posVect))
  full.readTable.posVect[empty.posVect] <- '*'
  
  if(length(readTable.list) == 2){
    full.readTable.posVect[names(readTable.list[[1]])] <- readTable.list[[1]]
    full.readTable.posVect[names(readTable.list[[2]])] <- readTable.list[[2]]
  }else{
    full.readTable.posVect[names(readTable.list[[1]])] <- readTable.list[[1]]
  }
  return(list(list(full.readTable.posVect)))
}

c4.phasing <- function(x, snpVect, snpDT){
  output.list <- list()
  
  snpNameVect <- names(snpVect)
  i <- 1
  for( current.snp in as.character(snpVect) ){
    
    if( current.snp %in% names(x) ){
      
      if( any( snpDT[,snpNameVect[i],with=F] == x[current.snp] ) ){
        output.list[[current.snp]] <- x[current.snp]
      }else if( x[current.snp] == '*' ){
        output.list[[current.snp]] <- '*'
      }else{
        output.list[[current.snp]] <- paste0('?',x[current.snp])
      }
      
    }else{
      output.list[[current.snp]] <- '-'
    }
    i <- i + 1
  }
  return(output.list)
}

c4.process_bamFile_to_samDT <- function( currentSmaple.id, rawData.dir ){
  cat('\n\t-- Processing',currentSample.id,'BAM file --')
  
  sample.dir <- list.files(path=file.path(rawData.dir,paste0('home/wmarin/C4Investigator/C4Investigator/output/tgp/')), pattern=currentSample.id,recursive = F,include.dirs = T,full.names = T)[[1]]
  sample.bam.path <- list.files(path=sample.dir,pattern = '.bam',recursive=T,full.names=T)
  
  cat('\n\t\tReading in',sample.bam.path)
  
  test.bam <- scanBam(sample.bam.path)
  currentSample.samDT <- as.data.table(test.bam[[1]])
  
  
  ## Use the reference name conversion list to apply common gene names
  currentSample.samDT$locus <- unlist(lapply(currentSample.samDT$rname, function(x){ referenceKeyList[[ as.character(x) ]] }))
  
  ## Pull out an index of reads that aligned to C4
  C4ReadsPosVect <- grep('C4',currentSample.samDT$locus,fixed=T)
  
  ## Initialize a column of the data table for storing whether the read aligned to C4 or not
  currentSample.samDT$isC4 <- F
  currentSample.samDT$isC4[C4ReadsPosVect] <- T ## Set the C4 bool to T for all C4 indexed reads
  
  currentSample.samDT$read_seq <- ''
  currentSample.samDT$read_seq <- sapply(currentSample.samDT$seq, function(x) as.character(x))
  currentSample.samDT[,seq:=NULL]
  currentSample.samDT[,qual:=NULL]
  
  
  ## Pull out the subset of C4 aligned reads as its own data table
  C4.currentSample.samDT <- copy( currentSample.samDT[isC4 == T][locus %in% c('C4A','C4B')] )
  rm(currentSample.samDT)
  
  ## Intialize and set a column for read lengths
  C4.currentSample.samDT$readLen <- nchar(C4.currentSample.samDT$read_seq)
  C4.currentSample.samDT[mrnm == 'hg38_knownGene_C4INSREG'][,pos := pos + 2860]
  
  cat('\n\t\tConverting alignment positions')
  ## Set the startPos and endPos for all reads
  C4.currentSample.samDT[,startPos := mapply(run.setStartPos, locus, pos)]
  C4.currentSample.samDT[,endPos := mapply(run.setEndPos, locus, pos, readLen)]
  
  ## Remove any positions that end up as NA (usually these are from reads that extend beyond the reference)
  C4.currentSample.samDT <- C4.currentSample.samDT[!is.na(C4.currentSample.samDT$endPos)]
  C4.currentSample.samDT <- C4.currentSample.samDT[!is.na(C4.currentSample.samDT$startPos)]
  
  ## Convert the SAM flags for all of the C4 aligned reads into an easily interpretable table
  bigSamFlagTable <- sapply(C4.currentSample.samDT$flag,samtable.flag_convert)
  
  C4.currentSample.samDT <- C4.currentSample.samDT[!unlist(bigSamFlagTable['notPrimaryAlignment',])]
  
  cat('\n\t\tConverting C4 aligned reads to read tables')
  C4.currentSample.samDT[, 'readTable' := alleleSetup.read_formatter( qname, read_seq, pos, locus, cigar), by=seq_len(nrow(C4.currentSample.samDT)) ]
  
  cat('\n\t\tDONE')
  return( C4.currentSample.samDT )
}

c4.samDT_to_phaseDT <- function( samDT ){
  cat('\n\t-- Preparing samDT for phasing --')
  cat('\n\t\tSelecting for reads between E25_1 - I28_20')
  toPhase.samDT <- copy( samDT[ startPos <= c4.feature.vect['I28_20'] & endPos >= c4.feature.vect['E25_1']] )
  
  cat('\n\t\tCombining paired-end read tables')
  toPhase.samDT$fullReadTable <- list()
  toPhase.samDT[, fullReadTable := c4.format_full_readTable(readTable),by=qname]
  toPhase.samDT$startPos <- as.integer(sapply(toPhase.samDT$fullReadTable, function(x) names(x)[1]))
  toPhase.samDT$endPos <- as.integer(sapply(toPhase.samDT$fullReadTable, function(x) names(x)[length(x)]))
  
  toPhase.samDT <- unique(toPhase.samDT, by=c('qname'))
  
  cat('\n\t\tDONE')
  return(toPhase.samDT)
}

c4.phaseDT_to_phaseMAT <- function(toPhase.samDT, toPhase.snpVect, toPhase.hetSnpDT){
  cat('\n\t-- Generating phase matrix --')
  cat('\n\t\tFinding reads spanning:',names(toPhase.snpVect))
  phase.mat <- t( sapply(toPhase.samDT$fullReadTable, c4.phasing, toPhase.snpVect, toPhase.hetSnpDT ) )
  colnames(phase.mat) <- names(toPhase.snpVect)
  
  phase.mat <- phase.mat[ apply( phase.mat, 1, function(x) sum(is_nuc(x)) > 1 ), ]
  cat('\n\t\tDONE')
  return(phase.mat)
}

AB.snpPosVect <- c('E26_129','E26_132','E26_140','E26_143','E26_145')
RgCh.snpPosVect <-  c('E28_111','E28_116','E28_125','E28_126')

c4.verify_AB_RgCh_phase <- function( phase.mat ){
  
  cat('\n\t-- Verifying phase of A/B and Rg/Ch --')
  
  C4A.snpVect <- c('C','G','T','G','C')
  C4B.snpVect <- c('T','C','A','C','T')
  names(C4A.snpVect) <- c('E26_129','E26_132','E26_140','E26_143','E26_145')
  names(C4B.snpVect) <- c('E26_129','E26_132','E26_140','E26_143','E26_145')
  
  C4Rg.snpVect <- c('G','T','T','C')
  C4Ch.snpVect <- c('C','C','G','G')
  names(C4Rg.snpVect) <- c('E28_111','E28_116','E28_125','E28_126')
  names(C4Ch.snpVect) <- c('E28_111','E28_116','E28_125','E28_126')
  
  phased.C4AB.posVect <- intersect(names(C4A.snpVect), colnames(phase.mat))
  phased.C4RgCh.posVect <- intersect(names(C4Rg.snpVect), colnames(phase.mat))
  
  C4AB.phase.mat <- unique( phase.mat[,phased.C4AB.posVect] )
  
  isAB.boolVect <- apply( C4AB.phase.mat, 1, function(x){
    isNuc.boolVect <- is_nuc( x )
    if( !any(isNuc.boolVect) ){return(T)}
    all( x[isNuc.boolVect] == C4A.snpVect[isNuc.boolVect] ) || all( x[isNuc.boolVect] == C4B.snpVect[isNuc.boolVect] )
  })
  if( all(isAB.boolVect) ){
    cat('\n\t\t- A/B verified in phase -')
  }else{
    cat('\n\t\t! A/B NOT in phase !')
    }
  
  C4RgCh.phase.mat <- unique( phase.mat[,phased.C4RgCh.posVect] )
  isRgCh.boolVect <- apply( C4RgCh.phase.mat, 1, function(x){
    isNuc.boolVect <- is_nuc( x )
    if( !any(isNuc.boolVect) ){return(T)}
    all( x[isNuc.boolVect] == C4Rg.snpVect[isNuc.boolVect] ) || all( x[isNuc.boolVect] == C4Ch.snpVect[isNuc.boolVect] )
  })
  if( all(isRgCh.boolVect) ){
    cat('\n\t\t- Rg/Ch verified in phase -')
  }else{
    cat('\n\t\t! Rg/Ch NOT in phase !')
    }
  
  cat('\n\t\tDONE')
  
  return( all(isRgCh.boolVect) & all(isAB.boolVect) )
}

c4.phase_AB_RgCh <- function( phase.mat ){
  cat('\n\t-- Phasing A/B with Rg/Ch --')
  
  toPhase.posVect <- c('E26_145','E28_111')
  
  temp.phaseMat <- t( apply(phase.mat[,toPhase.posVect], 1, gsub,pattern='*',replacement='_',fixed=T) )
  colnames(temp.phaseMat) <- toPhase.posVect
  
  temp.phaseMat <- temp.phaseMat[ apply( temp.phaseMat, 1, function(x) sum(is_nuc(x)) > 1 ), ]
  
  phase.convList <- list('E26_145'=list('C'='A','T'='B'),'E28_111'=list('G'='Rg','C'='Ch'))
  
  temp.phaseMat[,'E26_145'] <- as.character( phase.convList[['E26_145']][temp.phaseMat[,'E26_145']] )
  temp.phaseMat[,'E28_111'] <- as.character( phase.convList[['E28_111']][temp.phaseMat[,'E28_111']] )
  cat('\n\t\tDONE')
  return( table( apply( temp.phaseMat, 1, paste0, collapse='-' ) ) )
}

