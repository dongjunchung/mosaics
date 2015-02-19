
# generic methods for "MosaicsPeak" class

setMethod(
    f="show",
    signature="MosaicsPeak",
    definition=function( object ) {  
            
      peakParam <- object@peakParam
      
      analysisType <- peakParam@analysisType
      signalModel <- peakParam@signalModel
      FDR <- peakParam@FDR
      maxgap <- peakParam@maxgap
      minsize <- peakParam@minsize
      thres <- peakParam@thres
      
      cat( "Summary: MOSAiCS peak calling (class: MosaicsPeak)\n" )
      cat( "--------------------------------------------------\n" )
      cat( "final model: " )
      switch( analysisType,
          "OS" = {
              cat( "one-sample analysis " )
          },
          "TS" = {
              cat( "two-sample analysis (with M & GC) " )
          },
          "IO" = {
              cat( "two-sample analysis (input only) " )
          }
      )
      switch( signalModel,
          "1S" = {
              cat( "with one signal component\n" )
          },
          "2S" = {
              cat( "with two signal components\n" )
          }
      )        
      cat( "setting: FDR = ", FDR, ", maxgap = ", maxgap,
          ", minsize = ", minsize, ", thres = ", thres, "\n", sep="" )
      if ( length( object@peakList$peakStart ) > 0 ) {
          cat( "# of peaks = ", nrow(object@peakList), "\n", sep="" )
          cat( "median peak width = ", median(object@peakList$peakSize), "\n", sep="" )
          cat( "empirical FDR = ", round(empFDR(object)*10000)/10000, "\n", sep="" )
      } else {
          # exception handling (no peak case)
                  
          cat( "(no peak detected)\n" )
      }
    
      # information about read-level data
      
      if ( object@tagLoaded ) {
        cat( "read-level data is loaded.\n" )
        
        # suammry for read-level data
    
        nFrag <- unlist( lapply( read(object), function(x) length(x$ChIP) ) )
        sumRead <- sum(nFrag)
        medNumRead <- median(nFrag)
        
        cat( "ChIP sample:\n" )
        cat( "\tSequencing depth: ",seqDepth(object)[1],"\n", sep="" )
        cat( "\tNumber of utilized reads: ",sumRead,"\n", sep="" )
        cat( "\tMedian number of reads in each peak: ",medNumRead,"\n", sep="" )
        
        if ( !is.na(seqDepth(object)[2]) ) {
    
          nFrag <- unlist( lapply( read(object), function(x) length(x$Input) ) )
          sumRead <- sum(nFrag)
          medNumRead <- median(nFrag)
            
          cat( "Matched control sample:\n" )
          cat( "\tSequencing depth: ",seqDepth(object)[2],"\n", sep="" )
          cat( "\tNumber of utilized reads: ",sumRead,"\n", sep="" )
          cat( "\tMedian number of reads in each peak: ",medNumRead,"\n", sep="" )
        }
      } else {
        cat( "read-level data is currently not loaded.\n" )
        cat( "[Note] Read-level data is needed to run findSummit(), adjustBoundary(), and filterPeak().\n" ) 
        cat( "[Note] Run extractReads() to load read-level data.\n" )
      }
      
      cat( "--------------------------------------------------\n" )
    }
)

setMethod(
    f="print",
    signature="MosaicsPeak",
    definition=function( x ) {    
        peakList <- x@peakList
        
        if ( nrow(peakList) > 0 ) {
            printForm <- x@peakList
            return(printForm)
        } else {
            # exception handling (no peak case)
            
            message( "Info: no peak detected. Nothing can be returned." )           
        }
    }
)

setMethod(
    f="export",
    signature="MosaicsPeak",
    definition=function( object, type=NA, filename=NA ) {
      # error treatment: check invalid type
      
      if ( is.na(type) )
      {
        message( "Info: 'type' is not specified by the user." )
        message( "Info: 'type' is specified as 'txt' instead." )
        type <- "txt"        
      }
      
      allType <- c("txt","gff","bed","narrowPeak","broadPeak")
      invalidType <- TRUE
      for ( i in 1:length(type) )
      {
        if ( length(which(!is.na(match(type[i],allType))))==0 )
        {
          invalidType <- FALSE
        }
      }
      if ( !invalidType )
      {
        message( "Info: 'type' incorrect." )
        message( "Info: 'type' is specified as 'txt' instead." )
        type <- "txt"
      }
      
      # error treatment: 'filename' not specified
      
      if ( is.na(filename) )
      {
        message( "Info: 'filename' is not specified by the user." )
        message( "Info: 'filename' is specified as 'peakList' instead." )
        filename <- paste("peakList.",type,sep="")
      }
      
      # narrowPeak is supported only when tag data is loaded
      
      if ( type == "narrowPeak" & object@tagLoaded == FALSE ) {
        stop( "narrowPeak file format is supported only when read-level data is provided. Please run extractReads() first!" )
      }
    
      # no scientific notation
    
      scipenOrg <- options()$scipen
      options( scipen=999 )
        
      # export peak lists
      
      peakList <- object@peakList
      
      if ( type == "narrowPeak" | type == "broadPeak" ) {
        if ( object@tagLoaded == TRUE ) {
          summitSignal <- sapply( object@tagData@coverage, function(x) {
            if ( !is.na(x$ChIP[1,1]) ) {
              return( max(x$ChIP[,2]) )
            } else {
              return(0)
            }
          } )
        }
      } else {
        stop( "No read-level data provided. Please run extractReads() first!" )
      }
      
      if ( nrow(object@peakList) > 0 ) {        
        message( "Info: exporting the peak list..." )
        switch( type,
          "gff" = {
            .exportGFF( peakList=peakList, filename=filename )
          },
          "bed" = {
            .exportBED( peakList=peakList, filename=filename )
          },
          "txt" = {
            .exportTXT( peakList=peakList, filename=filename )
          },
          "narrowPeak" = {
            .exportEncodePeak( peakList=peakList, filename=filename, 
              fileformat="narrowPeak", summitSignal=summitSignal,
              inputExist=!is.na(seqDepth(object)[2]) )
          },
          "broadPeak" = {
            .exportEncodePeak( peakList=peakList, filename=filename, 
              fileformat="broadPeak", summitSignal=summitSignal,
              inputExist=!is.na(seqDepth(object)[2]) )
          }
        )
      } else {
        # exception handling (no peak case)
        
        message( "Info: no peak identifed. Nothing exported." )
      }
    
      # recover original setting for scientific notation
    
      options( scipen=scipenOrg )
  }
)

setMethod(
    f="plot",
    signature=c("MosaicsPeak","missing"),
    definition=function( x, y, filename=NA, peakNum=NA, ... ) {
      # supported only when tag data is loaded
      
      if ( x@tagLoaded == FALSE ) {
        stop( "No read-level data provided. Please run extractReads() first!" )
      }
      
      # generate PDF, if filename is provided
      
      if ( !is.na(filename) ) {
        pdf(filename)
      }
      
      # plot subset, if index is specified
      
      peakList <- print(x)
      nPeak <- nrow(peakList)
      
      if ( is.na(peakNum[1]) ) {
        peakNum <- 1:nPeak
      }
      
      # profile plots
      
      for ( j in peakNum ) {
        if ( !is.na( x@tagData@coverage[[j]]$ChIP[ 1, 1 ] ) ) {
          
          # ChIP track
          
          plot( x@tagData@coverage[[j]]$ChIP[,1], x@tagData@coverage[[j]]$ChIP[,2], 
            type="l", lwd=2, col="black",
            xlab="Genomic coordinates", ylab="Read count",
            main=names(x@tagData@coverage)[j] )
          
          # input track
          
          if ( !is.na( x@tagData@coverage[[j]]$Input[ 1, 1 ] ) ) {
            normC <- x@tagData@seqDepth[1] / x@tagData@seqDepth[2]
            
            lines( x@tagData@coverage[[j]]$Input[ , 1 ], 
              x@tagData@coverage[[j]]$Input[ , 2 ] * normC, 
              col = "gray", lwd = 2 )
            lines( x@tagData@coverage[[j]]$Input[ , 1 ], 
              x@tagData@coverage[[j]]$Input[, 2], 
              col = "pink", lwd = 2 )
            
            if ( colnames(peakList)[ ncol(peakList) ] == "summit" ) {
              legend( "topright", lwd=rep( 2, 5 ), lty=c( 1, 1, 1, 2, 2 ),
                col=c( "black", "gray", "pink", "black", "red" ),
                c( "ChIP data", "Control data (scaled)", "Control data (unscaled)", "Peak boundary", "Peak summit" ) )
            } else {
              legend( "topright", lwd=rep( 2, 4 ), lty=c( 1, 1, 1, 2 ),
                col=c( "black", "gray", "pink", "black" ),
                c( "ChIP data", "Control data (scaled)", "Control data (unscaled)", "Peak boundary" ) )
            }
          } else{
            if ( colnames(peakList)[ ncol(peakList) ] == "summit" ) {
              legend( "topright", lwd=rep( 2, 3 ), lty=c( 1, 2, 2 ),
                col=c( "black", "black", "red" ),
                c( "ChIP data", "Peak boundary", "Peak summit" ) )
            } else {
              legend( "topright", lwd=rep( 2, 2 ), lty=c( 1, 2 ),
                col=c( "black", "black" ),
                c( "ChIP data", "Peak boundary" ) )
            }
          }
          
          # peak boundaries
          
          abline( v = peakList[ j, 2 ], lty = 2, lwd = 2 )
          abline( v = peakList[ j, 3 ], lty = 2, lwd = 2 )
          
          # peak summits, if provided
          
          if ( colnames(peakList)[ ncol(peakList) ] == "summit" ) {
            abline( v = peakList[ j, ncol(peakList) ], col="red", lty = 2, lwd = 2 )
          }
        } else {
          # if there is no read in the region
          
          warning( "Peak plot cannot be provided for ",peakList[j,1],":",peakList[j,2],"-",peakList[j,3]," because there is no read in this peak region." )
        }
      }
      
      if ( !is.na(filename) ) {
        dev.off()
      }
    }
)

setMethod(
    f="empFDR",
    signature="MosaicsPeak",
    definition=function( object ) {
        return(object@empFDR)
    }
)

setMethod(
    f="bdBin",
    signature="MosaicsPeak",
    definition=function( object ) {
        return(object@bdBin)
    }
)

setMethod(
    f="coverage",
    signature="MosaicsPeak",
    definition=function( object ) {
      if ( object@tagLoaded == TRUE ) {  
        return(object@tagData@coverage)
      } else {
        stop( "Please run extractReads() first!" )
      }
    }
)

setMethod(
    f="read",
    signature="MosaicsPeak",
    definition=function( object ) {
      if ( object@tagLoaded == TRUE ) {  
        return(object@tagData@read)
      } else {
        stop( "Please run extractReads() first!" )
      }
    }
)

setMethod(
    f="seqDepth",
    signature="MosaicsPeak",
    definition=function( object ) {
      if ( object@tagLoaded == TRUE ) {  
        return(object@tagData@seqDepth)
      } else {
        stop( "Please run extractReads() first!" )
      }
    }
)
