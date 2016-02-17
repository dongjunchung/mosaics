
## Filter peaks from findSummit_trim_extend based on length and signal
## - minSignal = 20: Set to minRead*2
## - minLen = 200: Set to nucleosome length

setMethod(
  f="filterPeak",
  signature="MosaicsPeak",
  definition=function( object, 
    minSummitQuantile=0.01, minRead=10, FC=50, minSignal = 20, minLen = 200, 
    normC=NA, parallel=FALSE, nCore=8 ) {
    
    if ( object@tagLoaded == FALSE ) {
      stop( "[Note] Read-level data is needed to run filterPeak(). Run extractReads() to load read-level data." )
    }
    
    peakList <- print(object)
  
    # if not provided, estimate normC as ratio of seq depth of ChIP over Input
  
    if ( is.na(normC) ) {
      if ( !is.na(seqDepth(object)[2]) ) {
        normC <- seqDepth(object)[1] / seqDepth(object)[2]
      } else {
        normC <- 1
      }
    }
    
    # identify summits
    
    summit <- peakList[ , ncol(peakList) ]
    
    if ( parallel == TRUE ) {    
      summitSignalOrg <- mclapply( 1:nrow(peakList), function(j) {
        
        coverageChIP <- object@tagData@coverage[[j]]$ChIP
        
        if ( !is.na( coverageChIP[[1]] ) ) {
          xvar <- (coverageChIP[[1]]):(coverageChIP[[2]])
          yvar <- inverse.rle(coverageChIP[[3]])
          summitSignalChip <- yvar[ match( summit[j], xvar ) ]
        } else {
          summitSignalChip <- 0
        }
        if ( !is.na(seqDepth(object)[2]) && !is.na( object@tagData@coverage[[j]]$Input[[1]] ) ) {
          
          coverageInput <- object@tagData@coverage[[j]]$Input
          
          xvar <- (coverageInput[[1]]):(coverageInput[[2]])
          yvar <- inverse.rle(coverageInput[[3]])
          summitSignalInput <- yvar[ match( summit[j], xvar ) ]
          if ( is.na(summitSignalInput) ) {
            summitSignalInput <- 0
          }
        } else {
          summitSignalInput <- 0
        }
        return( c( summitSignalChip, summitSignalInput ) )
      }, mc.cores=nCore )
    } else {    
      summitSignalOrg <- lapply( 1:nrow(peakList), function(j) {
        
        coverageChIP <- object@tagData@coverage[[j]]$ChIP
        
        if ( !is.na( object@tagData@coverage[[j]]$ChIP[[1]] ) ) {
          xvar <- (coverageChIP[[1]]):(coverageChIP[[2]])
          yvar <- inverse.rle(coverageChIP[[3]])
          summitSignalChip <- yvar[ match( summit[j], xvar ) ]
        } else {
          summitSignalChip <- 0
        }
        if ( !is.na(seqDepth(object)[2]) && !is.na( object@tagData@coverage[[j]]$Input[[1]] ) ) {
          
          coverageInput <- object@tagData@coverage[[j]]$Input
          
          xvar <- (coverageInput[[1]]):(coverageInput[[2]])
          yvar <- inverse.rle(coverageInput[[3]])
          summitSignalInput <- yvar[ match( summit[j], xvar ) ]
          if ( is.na(summitSignalInput) ) {
            summitSignalInput <- 0
          }
        } else {
          summitSignalInput <- 0
        }
        return( c( summitSignalChip, summitSignalInput ) )
      } )
    }
    
    summitSignalChip <- sapply( summitSignalOrg, function(x) x[1] )
    summitSignalInput <- sapply( summitSignalOrg, function(x) x[2] )
    
    # determine summitCut
    
    summitCut <- max( quantile( summitSignalChip, minSummitQuantile ), 10 ) 
    
    # peak filtering #1
    
    improveCI <- 100 * ( summitSignalChip - summitSignalInput * normC ) / summitSignalChip
    improveCI[ is.na(improveCI) ] <- 0
      #NAs only occur when chip is 0.
    
    #numRead <- sapply( read(object), function(x) length(x$ChIP) )
    numRead <- object@tagData@numReads[,1]
    
    indRetained1 <- which( summitSignalChip >= summitCut & improveCI >= FC & numRead >= minRead )
    
    # peak filtering #2
    
    peakSize <- peakList[,3] - peakList[,2] + 1
    
    indRetained2 <- which( peakSize > minLen & summitSignalChip >= minSignal )    
    
    # filter peaks
    
    indRetained <- intersect( indRetained1, indRetained2 )
    
    object@peakList <- peakList[ indRetained, ]
    object@tagData@coverage <- object@tagData@coverage[ indRetained ]
    if ( object@tagData@keepReads == TRUE ) {
      object@tagData@read <- object@tagData@read[ indRetained ]
    }
    
    # summary
    
    cat( "------------------------------------------------------------\n" )
    cat( "Info: peak filtering summary\n" )
    cat( "------------------------------------------------------------\n" )
    cat( "# peaks before filtering: ", nrow(peakList), "\n" )
    cat( "# peaks after filtering #1: ", length(indRetained1), "\n" )
    cat( "# peaks after filtering #2: ", length(indRetained), "\n" )
    cat( "summitCut: ", summitCut, "\n" )
    cat( "------------------------------------------------------------\n" )
    
    return(object)
  }
)

