
## Trim and extend peak boundaries based on strengths of singal on boundaries
## summitSignal needs to be updated as boundaries are adjusted
## - FCoptimal <- 50
## - FCrelaxed <- 10

setMethod(
  f="adjustBoundary",
  signature="MosaicsPeak",
  definition=function( object, 
    minRead=10, extendFromSummit=100,
    trimMinRead1=1.5, trimFC1=5, extendMinRead1=2, extendFC1=50, 
    trimMinRead2=1.5, trimFC2=50, extendMinRead2=1.5, extendFC2=50,
    normC=NA, parallel=FALSE, nCore=8 ) {
    
      if ( object@tagLoaded == FALSE ) {
        stop( "[Note] Read-level data is needed to run adjustBoundary(). Run extractReads() to load read-level data." )
      }
    
      peakList <- print(object)
      chrCommon <- unique(peakList[,1])
      chrCommon <- sort(chrCommon)
    
      # if not provided, estimate normC as ratio of seq depth of ChIP over Input
    
      if ( is.na(normC) ) {
        if ( !is.na(seqDepth(object)[2]) ) {
          normC <- seqDepth(object)[1] / seqDepth(object)[2]
        } else {
          normC <- 1
        }
      }
    
      # adjust peak boundaries and update peak summits
    
      chipExist <- sapply( read(object), function(x) length(x$ChIP) ) > 0
      inputProvided <- !is.na(seqDepth(object)[2])
      if ( inputProvided ) {
        inputExist <- sapply( read(object), function(x) length(x$Input) ) > 0
      } else {
        inputExist <- rep( FALSE, length(read(object)) )
      }
      
      if ( parallel == TRUE ) {        
        out <- mclapply( 1:nrow(peakList), 
          function(j) .trimExtend( 
            chipExist=chipExist[j], inputProvided=inputProvided, inputExist=inputExist[j],
            chrID=peakList[ j, 1 ], peakStart=peakList[ j, 2 ], peakEnd=peakList[ j, 3 ], 
            summit=peakList[ j, ncol(peakList) ], stackedFragment=mosaics::coverage(object)[[j]], 
            normC=normC, extendFromSummit=extendFromSummit, minRead=minRead,
            trimMinRead1=trimMinRead1, trimFC1=trimFC1, extendMinRead1=extendMinRead1, extendFC1=extendFC1,
            trimMinRead2=trimMinRead2, trimFC2=trimFC2, extendMinRead2=extendMinRead2, extendFC2=extendFC2 ),
          mc.cores = nCore )
      } else {        
        out <- lapply( 1:nrow(peakList), 
          function(j) .trimExtend( 
            chipExist=chipExist[j], inputProvided=inputProvided, inputExist=inputExist[j],
            chrID=peakList[ j, 1 ], peakStart=peakList[ j, 2 ], peakEnd=peakList[ j, 3 ], 
            summit=peakList[ j, ncol(peakList) ], stackedFragment=mosaics::coverage(object)[[j]], 
            normC=normC, extendFromSummit=extendFromSummit, minRead=minRead,
            trimMinRead1=trimMinRead1, trimFC1=trimFC1, extendMinRead1=extendMinRead1, extendFC1=extendFC1,
            trimMinRead2=trimMinRead2, trimFC2=trimFC2, extendMinRead2=extendMinRead2, extendFC2=extendFC2 )
          )
      }
    
      # update original peak lists
    
      #peakList <- cbind( peakList, peakList[ , ncol(peakList) ] )
      #colnames(peakList)[ ( ncol(peakList) - 1 ) ] <- "summitOrg"
      #colnames(peakList)[ ncol(peakList) ] <- "summitAfterAdj"
    
      peakList[,2] <- sapply( out, function(x) x[1] )
      peakList[,3] <- sapply( out, function(x) x[2] )
      peakList[,4] <- peakList[,3] - peakList[,2] + 1
      peakList[,ncol(peakList)] <- sapply( out, function(x) x[3] )
    
      # summary
      
      cat( "------------------------------------------------------------\n" )
      cat( "Info: peak boundary adjustment summary\n" )
      cat( "------------------------------------------------------------\n" )
      cat( "# all peaks: ", nrow(object@peakList), "\n" )
      cat( "# peaks with trimmed boundaries: ", length(which( peakList[,2] > object@peakList[,2] | peakList[,3] < object@peakList[,3] )), "\n" )
      cat( "# peaks with extended boundaries: ", length(which( peakList[,2] < object@peakList[,2] | peakList[,3] > object@peakList[,3] )), "\n" )
      cat( "# peaks of which summits changed: ", length(which( peakList[,ncol(peakList)] != object@peakList[,ncol(object@peakList)] )), "\n" )
      cat( "------------------------------------------------------------\n" )
    
      object@peakList <- peakList
      return(object)
    }
)
