
# export all possible info

.exportTXT <- function( peakList, filename )
{
    outFormat <- peakList
    
    # Note (ver 1.5.5): chromStart 1-base & chromEnd inclusive.
    outFormat[ outFormat[,2] == 0, 2 ] <- 1
    	# first base should be 1, not 0, if we use 1-base system
        
    # variable names
    
    cat( file=filename )
    cat( as.character(colnames(outFormat)), file=filename, sep="\t", append=TRUE )
    cat( "\n", file=filename, append=TRUE )
    
    # peak list
     
    #for ( i in 1:nrow(outFormat) )
    #{
    #    cat( as.character(outFormat[i,]), file=filename, sep="\t", append=TRUE )
    #    cat( "\n", file=filename, append=TRUE )
    #}
    write.table( outFormat, file=filename, append=TRUE,
      sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE )
    
    message( "Info: peak file was exported in TXT format:" )
    message( "Info: file name = ", filename )
}

# GFF format (score = peak count)

.exportGFF <- function( peakList, filename )
{
    # GFF: seqname, source, feature, start, end, score, strand, frame, group
    # Note (ver 1.5.5): chromStart 1-base & chromEnd inclusive.
    
    outFormat <- data.frame( peakList$chrID, "MOSAiCS", "MOSAiCS_peak",
        peakList$peakStart, peakList$peakStop, peakList$aveChipCount,
        ".", ".", ".", stringsAsFactors=FALSE )
    outFormat[ outFormat[,4] <= 0, 4 ] <- 1
    	# first base should be 1, not 0, if we use 1-base system
    
    line0 <- 'track name=mosaicsPeaks description=\"MOSAiCS peaks\"'
    cat( as.character(line0), "\n", file=filename )
    #for ( i in 1:nrow(outFormat) )
    #{
    #    cat( as.character(outFormat[i,]), file=filename, sep="\t", append=TRUE )
    #    cat( "\n", file=filename, append=TRUE )
    #}
    write.table( outFormat, file=filename, append=TRUE,
      sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE )
    
    message( "Info: peak file was exported in GFF format:" )
    message( "Info: file name = ", filename )
}

# BED format (score = peak count)

.exportBED <- function( peakList, filename )
{
    # BED: (required) chrom, chromStart, chromEnd
    # BED: (optional) name, score, strand, thickStart, thinkEnd, itemRgb,
    #                   blockCount, blockSizes, blockStarts
    # Note (ver 1.5.5): chromStart 0-base & chromEnd NOT inclusive (i.e., 1-base).
    
    #outFormat <- data.frame( peakList$chrID,
    #    peakList$peakStart, peakList$peakStop, "MOSAiCS_peak", 
    #    peakList$aveChipCount, stringsAsFactors=FALSE )
    
    outFormat <- data.frame( peakList$chrID,
        (peakList$peakStart-1), peakList$peakStop, "MOSAiCS_peak", 
        peakList$aveChipCount, stringsAsFactors=FALSE )
    outFormat[ outFormat[,2] < 0, 2 ] <- 0
    	# first base becomes -1 by (peakList$peakStart-1)
        
    line0 <- 'track name=mosaicsPeaks description=\"MOSAiCS peaks\" useScore=1'
    cat( as.character(line0), "\n", file=filename )
    #for ( i in 1:nrow(outFormat) )
    #{
    #    cat( as.character(outFormat[i,]), file=filename, sep="\t", append=TRUE )
    #    cat( "\n", file=filename, append=TRUE )
    #}
    write.table( outFormat, file=filename, append=TRUE,
      sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE )
    
    message( "Info: peak file was exported in BED format:" )
    message( "Info: file name = ", filename )
}

# narrowPeak format

.exportEncodePeak <- function( peakList, filename, fileformat="broadPeak", summitSignal, inputExist )
{
    # narrowPeak: chrom, chromStart, chromEnd, name, score, strand, signalValue, pValue, qValue, peak
    
    #outFormat <- data.frame( peakList$chrID,
    #    peakList$peakStart, peakList$peakStop, "MOSAiCS_peak", 
    #    peakList$aveChipCount, stringsAsFactors=FALSE )
  
    if ( fileformat == "narrowPeak" ) {
      if ( inputExist ) {     
        # if control sample is provided, report ave2LogRatio
    
        outFormat <- data.frame( peakList$chrID,
          ( peakList$peakStart - 1 ), peakList$peakStop, "MOSAiCS_peak", 
          #peakList$aveLog2Ratio, ".", summitSignal, -log10(peakList$minP), -1, 
          #peakList$aveLog2Ratio, ".", summitSignal, peakList$logMinP, -1, 
          #peakList$aveLog2Ratio, ".", summitSignal, peakList$aveLogP, -1, 
          peakList$aveLog2Ratio, ".", summitSignal, peakList$logMinP, peakList$aveLogP,
          ( peakList$summit - peakList$peakStart ),
          stringsAsFactors=FALSE )
      } else {
        # if control sample is NOT provided, report aveChipCount
        
        outFormat <- data.frame( peakList$chrID,
          ( peakList$peakStart - 1 ), peakList$peakStop, "MOSAiCS_peak", 
          #peakList$aveChipCount, ".", summitSignal, -log10(peakList$minP), -1, 
          #peakList$aveChipCount, ".", summitSignal, peakList$logMinP, -1, 
          #peakList$aveChipCount, ".", summitSignal, peakList$aveLogP, -1, 
          peakList$aveChipCount, ".", summitSignal, peakList$logMinP, peakList$aveLogP,
          ( peakList$summit - peakList$peakStart ),
          stringsAsFactors=FALSE )      
      }
    } else if ( fileformat == "broadPeak" ) {
      if ( inputExist ) {     
        # if control sample is provided, report ave2LogRatio
    
        outFormat <- data.frame( peakList$chrID,
          ( peakList$peakStart - 1 ), peakList$peakStop, "MOSAiCS_peak", 
          #peakList$aveLog2Ratio, ".", summitSignal, -log10(peakList$minP), -1, 
          #peakList$aveLog2Ratio, ".", summitSignal, peakList$logMinP, -1, 
          #peakList$aveLog2Ratio, ".", summitSignal, peakList$aveLogP, -1, 
          peakList$aveLog2Ratio, ".", summitSignal, peakList$logMinP, peakList$aveLogP,
          stringsAsFactors=FALSE )
      } else {
        # if control sample is NOT provided, report aveChipCount
        
        outFormat <- data.frame( peakList$chrID,
          ( peakList$peakStart - 1 ), peakList$peakStop, "MOSAiCS_peak", 
          #peakList$aveChipCount, ".", summitSignal, -log10(peakList$minP), -1, 
          #peakList$aveChipCount, ".", summitSignal, peakList$logMinP, -1, 
          #peakList$aveChipCount, ".", summitSignal, peakList$aveLogP, -1, 
          peakList$aveChipCount, ".", summitSignal, peakList$logMinP, peakList$aveLogP,
          stringsAsFactors=FALSE )      
      }
    }
    outFormat[ outFormat[,2] < 0, 2 ] <- 0
      # first base becomes -1 by (peakList$peakStart-1)
        
    #line0 <- 'track name=mosaicsPeaks description=\"MOSAiCS peaks\" useScore=1'
    #cat( as.character(line0), "\n", file=filename )
    #for ( i in 1:nrow(outFormat) )
    #{
    #    cat( as.character(outFormat[i,]), file=filename, sep="\t", append=TRUE )
    #    cat( "\n", file=filename, append=TRUE )
    #}
    write.table( outFormat, file=filename, append=FALSE,
      sep="\t", quote=FALSE, row.names=FALSE, col.names=FALSE )
    
    if ( fileformat == "narrowPeak" ) {
      message( "Info: peak file was exported in narrowPeak format:" )
    } else if ( fileformat == "broadPeak" ) {
      message( "Info: peak file was exported in broadPeak format:" )
    }
    message( "Info: file name = ", filename )
}
