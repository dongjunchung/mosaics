
## Extract reads corresponding to each peak region

setMethod(
  f="extractReads",
  signature="MosaicsPeak",
  definition=function( object, chipFile=NULL, chipFileFormat=NULL,
    chipPET=FALSE, chipFragLen=200,
    controlFile=NULL, controlFileFormat=NULL, 
    controlPET=FALSE, controlFragLen=200,
    parallel=FALSE, nCore=8, tempDir=NULL, perl="perl" )
  { 
    
    keepReads=TRUE
      # keep read-level data
    
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
    
    # rearranage results
    
    if ( !is.null(controlFile) ) {
      seqDepth <- c( outChIP$seqDepth, outInput$seqDepth )
    } else {
      seqDepth <- c( outChIP$seqDepth, NA )
    }
    
    fragSet <- stackedFragment <- vector( "list", nPeak )
    
    for ( i in 1:nPeak ) {
      fragSet[[i]] <- vector( "list", 2 )
      stackedFragment[[i]] <- vector( "list", 2 )  
      
      if( !is.null(outChIP$fragSet[[i]]) ) {  
        fragSet[[i]]$ChIP <- outChIP$fragSet[[i]]
        stackedFragment[[i]]$ChIP <- outChIP$stackedFragment[[i]]
      } else {
        fragSet[[i]]$ChIP <- GRanges()
        stackedFragment[[i]]$ChIP <- matrix(NA)
      }
      
      if ( !is.null(controlFile) && !is.null(outInput$fragSet[[i]]) ) {
        fragSet[[i]]$Input <- outInput$fragSet[[i]]
        stackedFragment[[i]]$Input <- outInput$stackedFragment[[i]]
      } else {
        fragSet[[i]]$Input <- GRanges()
        stackedFragment[[i]]$Input <- matrix(NA)
      }
    }
    
    names(fragSet) <- names(stackedFragment) <- 
      paste( peakList[,1], ":", peakList[,2], "-", peakList[,3], sep="" )
    
    # info about preprocessing
    
    cat( "------------------------------------------------------------\n" )
    cat( "Info: Preprocessing summary\n" )
    cat( "------------------------------------------------------------\n" )
    cat( "Number of chromosomes: ",length(chrCommon),"\n", sep="" )
    cat( "Number of peaks: ",nPeak,"\n", sep="" )
    
    nFrag <- unlist( lapply( fragSet, function(x) length(x$ChIP) ) )
    sumRead <- sum(nFrag)
    medNumRead <- median(nFrag)
    
    cat( "ChIP sample:\n" )
    cat( "\tTag type: ",ifelse(chipPET,"PET","SET"),"\n", sep="" )
    cat( "\tSequencing depth: ",seqDepth[1],"\n", sep="" )
    cat( "\tNumber of utilized reads: ",sumRead,"\n", sep="" )
    cat( "\tMedian number of reads in each peak: ",medNumRead,"\n", sep="" )
    if ( !is.null(controlFile) ) {
      
      nFrag <- unlist( lapply( fragSet, function(x) length(x$Input) ) )
      sumRead <- sum(nFrag)
      medNumRead <- median(nFrag)
        
      cat( "Matched control sample:\n" )
      cat( "\tTag type: ",ifelse(controlPET,"PET","SET"),"\n", sep="" )
      cat( "\tSequencing depth: ",seqDepth[2],"\n", sep="" )
      cat( "\tNumber of utilized reads: ",sumRead,"\n", sep="" )
      cat( "\tMedian number of reads in each peak: ",medNumRead,"\n", sep="" )
    }
    
    cat( "------------------------------------------------------------\n" )
    
    # update object
    
    object@tagLoaded <- TRUE
    object@tagData <- new( "TagData", 
      read=fragSet, coverage=stackedFragment, seqDepth=seqDepth )    
    return(object)
  }
)
