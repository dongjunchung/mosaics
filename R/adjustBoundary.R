
## Trim and extend peak boundaries based on strengths of singal on boundaries
## summitSignal needs to be updated as boundaries are adjusted
## - FCoptimal <- 50
## - FCrelaxed <- 10

setMethod(
  f="adjustBoundary",
  signature="MosaicsPeak",
  definition=function( object, minRead=10, extendFromSummit=100,
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
    
      #chipExist <- sapply( read(object), function(x) length(x$ChIP) ) > 0
      chipExist <- object@tagData@numReads[,1] > 0
      inputProvided <- !is.na(seqDepth(object)[2])
      if ( inputProvided ) {
        #inputExist <- sapply( read(object), function(x) length(x$Input) ) > 0
        inputExist <- object@tagData@numReads[,2] > 0
      } else {
        #inputExist <- rep( FALSE, length(read(object)) )
        inputExist <- rep( FALSE, nrow(object@tagData@numReads) )
      }
      
      if ( parallel == TRUE ) {        
        out <- mclapply( 1:nrow(peakList), 
          function(j) .trimExtend( 
            chipExist=chipExist[j], inputProvided=inputProvided, inputExist=inputExist[j],
            chrID=peakList[ j, 1 ], peakStart=peakList[ j, 2 ], peakEnd=peakList[ j, 3 ], 
            summitSignal=peakList[ j, (ncol(peakList)-1) ], summit=peakList[ j, ncol(peakList) ], 
            stackedFragment=object@tagData@coverage[[j]], 
            normC=normC, extendFromSummit=extendFromSummit, minRead=minRead,
            trimMinRead1=trimMinRead1, trimFC1=trimFC1, extendMinRead1=extendMinRead1, extendFC1=extendFC1,
            trimMinRead2=trimMinRead2, trimFC2=trimFC2, extendMinRead2=extendMinRead2, extendFC2=extendFC2 ),
          mc.cores = nCore )
      } else {        
        out <- lapply( 1:nrow(peakList), 
          function(j) .trimExtend( 
            chipExist=chipExist[j], inputProvided=inputProvided, inputExist=inputExist[j],
            chrID=peakList[ j, 1 ], peakStart=peakList[ j, 2 ], peakEnd=peakList[ j, 3 ], 
            summitSignal=peakList[ j, (ncol(peakList)-1) ], summit=peakList[ j, ncol(peakList) ], 
            stackedFragment=object@tagData@coverage[[j]], 
            normC=normC, extendFromSummit=extendFromSummit, minRead=minRead,
            trimMinRead1=trimMinRead1, trimFC1=trimFC1, extendMinRead1=extendMinRead1, extendFC1=extendFC1,
            trimMinRead2=trimMinRead2, trimFC2=trimFC2, extendMinRead2=extendMinRead2, extendFC2=extendFC2 )
          )
      }
    
      # update peak annotations

      peakStart <- sapply( out, function(x) x[1] )
      peakStop <- sapply( out, function(x) x[2] )
    
      loc_list <- split( 1:nrow(peakList), peakList[,1] )
      peakStart_list <- split( peakStart, peakList[,1] )
      peakStop_list <- split( peakStop, peakList[,1] )
    
      betapH_list <- split( object@postProb, object@chrID )
      coord_list <- split( object@coord, object@chrID )
      Y_list <- split( object@tagCount, object@chrID )
      switch( object@peakParam@analysisType,
        OS = {
          M_list <- split( object@mappability, object@chrID )
          GC_list <- split( object@gcContent, object@chrID )
          nRatio <- 1
        },
        TS = {
          X_list <- split( object@input, object@chrID )
          M_list <- split( object@mappability, object@chrID )
          GC_list <- split( object@gcContent, object@chrID )
          nRatio <- object@seqDepth[1] / object@seqDepth[2]
        },
        IO = {
          X_list <- split( object@input, object@chrID )
          nRatio <- object@seqDepth[1] / object@seqDepth[2]
        }
      ) 
    
      #chrList <- sort(unique(object@chrID))
      chrList <- as.character(chrCommon)
      
      for ( chr in 1:length(chrList) ) {
        # extract data for given chromosome
        
        loc_chr <- loc_list[[ chrList[chr] ]]
        
        peakStart_chr <- peakStart_list[[ chrList[chr] ]]
        peakStop_chr <- peakStop_list[[ chrList[chr] ]]
        
        betapH_chr <- betapH_list[[ chrList[chr] ]]
        coord_chr <- coord_list[[ chrList[chr] ]]
        Y_chr <- Y_list[[ chrList[chr] ]]
        switch( object@peakParam@analysisType,
            OS = {
                X_chr <- NA
                M_chr <- M_list[[ chrList[chr] ]]
                GC_chr <- GC_list[[ chrList[chr] ]]
            },
            TS = {
                X_chr <- X_list[[ chrList[chr] ]]
                M_chr <- M_list[[ chrList[chr] ]]
                GC_chr <- GC_list[[ chrList[chr] ]]
            },
            IO = {
                X_chr <- X_list[[ chrList[chr] ]]
                M_chr <- NA
                GC_chr <- NA
            }
        ) 
                
        final_peakset_chr <- .annotatePeak( 
          peakStart_chr=peakStart_chr, peakStop_chr=peakStop_chr, 
          coord_chr=coord_chr, analysisType=object@peakParam@analysisType,
          Y_chr=Y_chr, X_chr=X_chr, M_chr=M_chr, GC_chr=GC_chr, pp_chr=betapH_chr[,3], 
          nRatio=nRatio )
        
        peakList[ loc_chr, 2:(ncol(peakList)-2) ] <- final_peakset_chr
      }
    
      # update peak regions
    
      #peakList[,2] <- peakStart
      #peakList[,3] <- peakStop
      #peakList[,4] <- peakList[,3] - peakList[,2] + 1
      peakList[,(ncol(peakList)-1)] <- sapply( out, function(x) x[3] )
      peakList[,ncol(peakList)] <- sapply( out, function(x) x[4] )
        
      # summary
      
      cat( "------------------------------------------------------------\n" )
      cat( "Info: peak boundary adjustment summary\n" )
      cat( "------------------------------------------------------------\n" )
      cat( "# all peaks: ", nrow(object@peakList), "\n" )
      cat( "# peaks with trimmed boundaries: ", length(which( peakList[,2] > object@peakList[,2] | peakList[,3] < object@peakList[,3] )), "\n" )
      cat( "# peaks with extended boundaries: ", length(which( peakList[,2] < object@peakList[,2] | peakList[,3] > object@peakList[,3] )), "\n" )
      cat( "# peaks of which summits changed: ", length(which( peakList[,ncol(peakList)] != object@peakList[,ncol(object@peakList)] )), "\n" )
      cat( "normC: ", normC, "\n" )
      cat( "------------------------------------------------------------\n" )
    
      object@peakList <- peakList
      return(object)
    }
)
