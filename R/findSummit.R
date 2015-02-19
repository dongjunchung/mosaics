
## Identify a summit in each peak region

setMethod(
  f="findSummit",
  signature="MosaicsPeak",
  definition=function( object, parallel=FALSE, nCore=8 ) {
    
    if ( object@tagLoaded == FALSE ) {
      stop( "[Note] Read-level data is needed to run findSummit(). Run extractReads() to load read-level data." )
    }
      
    # identify a summit for each peak region
    # [Note] summit should exist within the peak region
    
    peakList <- print(object)
    peakStart <- peakList[,2]
    peakEnd <- peakList[,3]
    
    if ( parallel == TRUE ) {        
      summits <- mclapply( as.list(1:nrow(peakList)), function(j) {
          
        if ( !is.na(coverage(object)[[j]]$ChIP[1,1]) ) {
          
          profileChip <- coverage(object)[[j]]$ChIP
          locStart <- match( peakStart[j], profileChip[,1] )
          locEnd <- match( peakEnd[j], profileChip[,1] )
          if ( is.na(locStart) ) {
            locStart <- 1
          }
          locEnd <- match( peakEnd[j], profileChip[,1] )
          if ( is.na(locEnd) ) {
            locEnd <- nrow(profileChip)
          }
        
          locSummit <- which( profileChip[ locStart:locEnd, 2 ] == max( profileChip[ locStart:locEnd, 2 ] ) )
          if ( length(locSummit) > 1 ) {
            locSummit <- round(mean(locSummit))
          }
          return( profileChip[ locStart:locEnd, 1 ][ locSummit ] )
          
        } else {
          # if there is no read, simply return midpoint of the peak region
          
          return( round( ( peakStart[j] + peakEnd[j] ) / 2 ) )
          
        }
        
      }, mc.cores = nCore )
    } else {        
      summits <- lapply( as.list(1:nrow(peakList)), function(j) {
          
        if ( !is.na(coverage(object)[[j]]$ChIP[1,1]) ) {
          
          profileChip <- coverage(object)[[j]]$ChIP
          locStart <- match( peakStart[j], profileChip[,1] )
          if ( is.na(locStart) ) {
            locStart <- 1
          }
          locEnd <- match( peakEnd[j], profileChip[,1] )
          if ( is.na(locEnd) ) {
            locEnd <- nrow(profileChip)
          }
        
          locSummit <- which( profileChip[ locStart:locEnd, 2 ] == max( profileChip[ locStart:locEnd, 2 ] ) )
          if ( length(locSummit) > 1 ) {
            locSummit <- round(mean(locSummit))
          }
          return( profileChip[ locStart:locEnd, 1 ][ locSummit ] )
          
        } else {
          # if there is no read, simply return midpoint of the peak region
          
          return( round( ( peakStart[j] + peakEnd[j] ) / 2 ) )
          
        }
        
      } )
    }
  
    # add summit information as the last column in the peak list
  
    object@peakList <- cbind( object@peakList, unlist(summits) )
    #object@peakList <- cbind( object@peakList, summits - object@peakList[,2] )
    #  # summits are represented as 0-based
    colnames(object@peakList)[ length(object@peakList) ] <- "summit"
  
    return(object)
  } 
)
