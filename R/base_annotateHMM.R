.annotateHMM <- function( object, analysisType,
	maxgap=0, minsize=0, thres=0, binsize=200, nRatio=1 ) {
		
    inputdata <- object$inputdata
    bd_bin <- object$bd_bin
    betapH <- object$betapH
		
	# extract data for given chromosome
    
    coord_chr <- inputdata[,1]
    peakcall_chr <- inputdata[,2]
    gMat_chr <- rbind( inputdata[,3], inputdata[,4] )
    Y_chr <- inputdata[,5]
    
    switch( analysisType,
        OS = {
            X_chr <- NA
            M_chr <- inputdata[,6]
            GC_chr <- inputdata[,7]
        },
        TS = {
            X_chr <- inputdata[,6]
            M_chr <- inputdata[,7]
            GC_chr <- inputdata[,8]
        },
        IO = {
            X_chr <-inputdata[,6]
            M_chr <- NA
            GC_chr <- NA
        }
    ) 
	
	  if ( length(which( bd_bin == 1 & Y_chr >= thres )) > 0 ) {
		  #################### if we still have peaks after tag count thresholding
		
		  bd_ID <- which( bd_bin == 1 & Y_chr >= thres )
    
	    # update empirical FDR
	    
	    empFDR <- sum( betapH[ bd_ID ] ) / length(bd_ID)
	    
	    # process peaks
	    
	    peak_range <- IRanges( start=coord_chr[ bd_ID ], 
	        end=coord_chr[ bd_ID ] + binsize - 1 )
	    peak_reduced <- reduce(peak_range)
	    peak_start <- start(peak_reduced)
		  peak_stop <- end(peak_reduced)
	    
	    # merge close peaks if distance<=maxgap
	    
	    coord_org <- IRanges( start=peak_start, end=peak_stop+maxgap )
	    coord_merged <- reduce(coord_org)
	    peak_start <- start(coord_merged)
	    peak_stop <- end(coord_merged) - maxgap                    
	    
	    # filter peaks smaller than minsize & order by coordinates
	    
	    peaksize <- peak_stop - peak_start + 1
	    filterID = which( peaksize <= minsize )
	    if ( length(filterID) > 0 )
	    {
	        peak_start <- peak_start[ -filterID ]
	        peak_stop <- peak_stop[ -filterID ]
	    }
	    peak_start <- peak_start[ order(peak_start) ]        
	    peak_stop <- peak_stop[ order(peak_start) ]
	    peaksize <- peak_stop - peak_start + 1
	    
	    # calculate additional info
      
      final_peakset_chr <- .annotatePeak( 
        peakStart_chr=peak_start, peakStop_chr=peak_stop, coord_chr=coord_chr, analysisType=analysisType,
        Y_chr=Y_chr, X_chr=X_chr, M_chr=M_chr, GC_chr=GC_chr, pp_chr=betapH, nRatio=nRatio )
    } else {
        #################### if there is no peak
        
        final_peakset_chr=data.frame()
        empFDR=0
    }
	
	return( list( final_peakset=final_peakset_chr, empFDR=empFDR ) )
}