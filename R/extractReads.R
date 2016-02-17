
## Extract reads corresponding to each peak region

setMethod(
  f="extractReads",
  signature="MosaicsPeak",
  definition=function( object, chipFile=NULL, chipFileFormat=NULL,
    chipPET=FALSE, chipFragLen=200,
    controlFile=NULL, controlFileFormat=NULL, 
    controlPET=FALSE, controlFragLen=200, keepReads=FALSE,
    parallel=FALSE, nCore=8, tempDir=NULL, perl="perl" )
  { 
    
    # summarize peak info
    
    peakList <- print(object)
    nPeak <- nrow(peakList)
    chrCommon <- unique(peakList[,1])
    chrCommon <- sort(chrCommon)
    
    # process read files
    
    message( "Info: Loading and processing ChIP sample file..." )
    
    outChIP <- .loadReadData( object=object, readfile=chipFile,
      fileFormat=chipFileFormat, PET=chipPET, fragLen=chipFragLen, keepReads=keepReads,
      parallel=parallel, nCore=nCore, tempDir=tempDir, perl=perl )    
    
    if ( !is.null(controlFile) ) {
      message( "Info: Loading and processing matched control sample file..." )
      
      outInput <- .loadReadData( object=object, readfile=controlFile,
        fileFormat=controlFileFormat, PET=controlPET, fragLen=controlFragLen, keepReads=keepReads,
        parallel=parallel, nCore=nCore, tempDir=tempDir, perl=perl )
    }
    
    # rearranage results: seqDepth
    
    if ( !is.null(controlFile) ) {
      seqDepth <- c( outChIP$seqDepth, outInput$seqDepth )
    } else {
      seqDepth <- c( outChIP$seqDepth, NA )
    }
    
    # rearranage results: stackedFragment
    
    stackedFragment <- vector( "list", nPeak )
    
    for ( i in 1:nPeak ) {
      stackedFragment[[i]] <- vector( "list", 2 )  
      
      if( !is.na(outChIP$stackedFragment[[i]][[1]]) ) {  
        stackedFragment[[i]]$ChIP <- outChIP$stackedFragment[[i]]
      } else {
        stackedFragment[[i]]$ChIP <- vector( "list", 3 )
        stackedFragment[[i]]$ChIP[[1]] <- stackedFragment[[i]]$ChIP[[2]] <- 
          stackedFragment[[i]]$ChIP[[3]] <- NA
      }
      
      if ( !is.null(controlFile) && !is.na(outInput$stackedFragment[[i]][[1]]) ) {
        stackedFragment[[i]]$Input <- outInput$stackedFragment[[i]]
      } else {
        stackedFragment[[i]]$Input <- vector( "list", 3 )
        stackedFragment[[i]]$Input[[1]] <- stackedFragment[[i]]$Input[[2]] <- 
          stackedFragment[[i]]$Input[[3]] <- NA
      }
    }
    
    names(stackedFragment) <- paste( peakList[,1], ":", peakList[,2], "-", peakList[,3], sep="" )
    
    # rearranage results: fragSet
    
    if ( keepReads == TRUE ) {
      fragSet <- vector( "list", nPeak )
      
      for ( i in 1:nPeak ) {
        fragSet[[i]] <- vector( "list", 2 )
        
        if( !is.na(outChIP$stackedFragment[[i]][[1]]) ) {  
          fragSet[[i]]$ChIP <- outChIP$fragSet[[i]]
        } else {
          fragSet[[i]]$ChIP <- GRanges()
        }
        
        if ( !is.null(controlFile) && !is.na(outInput$stackedFragment[[i]][[1]]) ) {
          fragSet[[i]]$Input <- outInput$fragSet[[i]]
        } else {
          fragSet[[i]]$Input <- GRanges()
        }
      }
      
      names(fragSet) <- paste( peakList[,1], ":", peakList[,2], "-", peakList[,3], sep="" )
    } else {
      fragSet <- list()
    }
    
    # rearranage results: numReads
    
    numReads <- matrix( NA, nPeak, 2 )
    numReads[,1] <- outChIP$numReads
    if ( !is.null(controlFile) ) {
      numReads[,2] <- outInput$numReads
    }
    rownames(numReads) <- paste( peakList[,1], ":", peakList[,2], "-", peakList[,3], sep="" )
    colnames(numReads) <- c( "ChIP", "Control" )
    
    # info about preprocessing
    
    cat( "------------------------------------------------------------\n" )
    cat( "Info: Preprocessing summary\n" )
    cat( "------------------------------------------------------------\n" )
    cat( "Number of chromosomes: ",length(chrCommon),"\n", sep="" )
    cat( "Number of peaks: ",nPeak,"\n", sep="" )
    sumRead <- sum(numReads[,1])
    medNumRead <- median(numReads[,1])
    
    cat( "ChIP sample:\n" )
    cat( "\tTag type: ",ifelse(chipPET,"PET","SET"),"\n", sep="" )
    cat( "\tSequencing depth: ",seqDepth[1],"\n", sep="" )
    cat( "\tNumber of utilized reads: ",sumRead,"\n", sep="" )
    cat( "\tMedian number of reads in each peak: ",medNumRead,"\n", sep="" )
    if ( !is.null(controlFile) ) {
      sumRead <- sum(numReads[,2])
      medNumRead <- median(numReads[,2])
        
      cat( "Matched control sample:\n" )
      cat( "\tTag type: ",ifelse(controlPET,"PET","SET"),"\n", sep="" )
      cat( "\tSequencing depth: ",seqDepth[2],"\n", sep="" )
      cat( "\tNumber of utilized reads: ",sumRead,"\n", sep="" )
      cat( "\tMedian number of reads in each peak: ",medNumRead,"\n", sep="" )
    }
    
    cat( "------------------------------------------------------------\n" )
    
    # update object
    
    object@tagLoaded <- TRUE
    #object@tagData <- new( "TagData", 
    #  read=fragSet, coverage=stackedFragment, seqDepth=seqDepth )    
    object@tagData <- new( "TagData", 
      read=fragSet, numReads=numReads, coverage=stackedFragment, keepReads=keepReads ) 
    object@seqDepth <- seqDepth
    return(object)
  }
)
