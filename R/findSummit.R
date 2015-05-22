
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
          
        if ( !is.na(object@tagData@coverage[[j]]$ChIP[[1]]) ) {
          
          xvar <- (object@tagData@coverage[[j]]$ChIP[[1]]):(object@tagData@coverage[[j]]$ChIP[[2]])
          yvar <- inverse.rle(object@tagData@coverage[[j]]$ChIP[[3]])
          
          profileChip <- cbind( xvar, yvar )
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
            firstPeak <- which.max(profileChip[ locStart:locEnd, 2 ])
            firstBreak <- locSummit[ which( diff(locSummit) != 1 )[1] ]
            if ( is.na(firstBreak) ) {
              firstBreak <- max(locSummit)
            }
            firstBlock <- firstPeak:firstBreak
            locSummit <- floor(mean(firstBlock))
          }
          
          summitCoord <- profileChip[ locStart:locEnd, 1 ][ locSummit ]
          summitSignal <- profileChip[ locStart:locEnd, 2 ][ locSummit ] 
          
          return( c( summitCoord, summitSignal ) )
          
        } else {
          # if there is no read, simply return midpoint of the peak region
          
          summitCoord <- floor( ( peakStart[j] + peakEnd[j] ) / 2 )
          summitSignal <- 0
          
          return( c( summitCoord, summitSignal ) )
        }
        
      }, mc.cores = nCore )
    } else {        
      summits <- lapply( as.list(1:nrow(peakList)), function(j) {
          
        if ( !is.na(object@tagData@coverage[[j]]$ChIP[[1]]) ) {
          
          xvar <- (object@tagData@coverage[[j]]$ChIP[[1]]):(object@tagData@coverage[[j]]$ChIP[[2]])
          yvar <- inverse.rle(object@tagData@coverage[[j]]$ChIP[[3]])
          
          profileChip <- cbind( xvar, yvar )
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
            firstPeak <- which.max(profileChip[ locStart:locEnd, 2 ])
            firstBreak <- locSummit[ which( diff(locSummit) != 1 )[1] ]
            if ( is.na(firstBreak) ) {
              firstBreak <- max(locSummit)
            }
            firstBlock <- firstPeak:firstBreak
            locSummit <- floor(mean(firstBlock))
          }
          
          summitCoord <- profileChip[ locStart:locEnd, 1 ][ locSummit ]
          summitSignal <- profileChip[ locStart:locEnd, 2 ][ locSummit ] 
          
          return( c( summitCoord, summitSignal ) )
          
        } else {
          # if there is no read, simply return midpoint of the peak region
          
          summitCoord <- floor( ( peakStart[j] + peakEnd[j] ) / 2 )
          summitSignal <- 0
          
          return( c( summitCoord, summitSignal ) )
          
        }
        
      } )
    }
  
    # add summit information as the last column in the peak list
    
    summitCoord <- sapply( summits, function(x) x[1] )
    summitSignal <- sapply( summits, function(x) x[2] )
  
    object@peakList <- cbind( object@peakList, summitSignal, summitCoord )
    #object@peakList <- cbind( object@peakList, summits - object@peakList[,2] )
    #  # summits are represented as 0-based
    
    colnames(object@peakList)[ (length(object@peakList)-1):length(object@peakList) ] <- c( "summitSignal", "summit" )
  
    return(object)
  } 
)
