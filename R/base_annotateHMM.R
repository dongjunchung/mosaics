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
	        
	    peak_range <- IRanges( start=peak_start, end=peak_stop )
	    info_range <- IRanges( start=coord_chr, end=coord_chr )
	    mm <- as.matrix( findOverlaps( peak_range, info_range ) )
	    matchlist <- split( mm[,2], mm[,1] )
	    
	    switch( analysisType,
	        OS = {        
	            # extract info
	            
	            peak_info <- sapply( matchlist, 
	                function(x) {
	                    betapH_i <- betapH[x]
	                    aveP <- mean(betapH_i)
	                    minP <- min(betapH_i)
	                    aveChipCount <- mean(Y_chr[x])
	                    maxChipCount <- max(Y_chr[x])
	                    aveM <- mean(M_chr[x])
	                    aveGC <- mean(GC_chr[x]) 
	                    return( c( aveP, minP, aveChipCount, maxChipCount, aveM, aveGC ) )
	                }
	            )
	                
	            peak_info <- t(peak_info)
	            
	            # combine all
	                    
	            final_peakset_chr <- data.frame(
	                peak_start, peak_stop, peaksize, peak_info, 
	                stringsAsFactors=FALSE )
	            colnames(final_peakset_chr) <-
	                c('peakStart','peakStop','peakSize','aveP','minP',
	                'aveChipCount','maxChipCount','map','GC')
	        },
	        TS = {        
	            # extract info
	        
	            peak_info <- sapply( matchlist, 
	                function(x) {
	                    betapH_i <- betapH[x]
	                    aveP <- mean(betapH_i)
	                    minP <- min(betapH_i)
	                    aveChipCount <- mean(Y_chr[x])
	                    maxChipCount <- max(Y_chr[x])
	                    aveInputCount <- mean(X_chr[x])
	                    aveInputCountScaled <- aveInputCount * nRatio
	                    aveLog2Ratio <- mean( log2( (Y_chr[x]+1) / (X_chr[x]*nRatio+1) ) )
	                    aveM <- mean(M_chr[x])
	                    aveGC <- mean(GC_chr[x]) 
	                    return( c( aveP, minP, aveChipCount, maxChipCount, 
	                        aveInputCount, aveInputCountScaled, aveLog2Ratio, aveM, aveGC ) )
	                }
	            )
	                
	            peak_info <- t(peak_info)
	            
	            # combine all
	                    
	            final_peakset_chr <- data.frame( 
	                peak_start, peak_stop, peaksize, peak_info, 
	                stringsAsFactors=FALSE )
	            colnames(final_peakset_chr) <-
	                c('peakStart','peakStop','peakSize','aveP','minP',
	                'aveChipCount','maxChipCount',
	                'aveInputCount','aveInputCountScaled','aveLog2Ratio','map','GC')
	        },
	        IO = {        
	            # extract info
	            
	            peak_info <- sapply( matchlist, 
	                function(x) {
	                    betapH_i <- betapH[x]
	                    aveP <- mean(betapH_i)
	                    minP <- min(betapH_i)
	                    aveChipCount <- mean(Y_chr[x])
	                    maxChipCount <- max(Y_chr[x])
	                    aveInputCount <- mean(X_chr[x])
	                    aveInputCountScaled <- aveInputCount * nRatio
	                    aveLog2Ratio <- mean( log2( (Y_chr[x]+1) / (X_chr[x]*nRatio+1) ) )
	                    return( c( aveP, minP, aveChipCount, maxChipCount, 
	                        aveInputCount, aveInputCountScaled, aveLog2Ratio ) )
	                }
	            )
	                
	            peak_info <- t(peak_info)
	        
	            # combine all
	                    
	            final_peakset_chr <- data.frame( 
	                peak_start, peak_stop, peaksize, peak_info, 
	                stringsAsFactors=FALSE )
	            colnames(final_peakset_chr) <-
	                c('peakStart','peakStop','peakSize','aveP','minP',
	                'aveChipCount','maxChipCount',
	                'aveInputCount','aveInputCountScaled','aveLog2Ratio')
	        }
	    )
    } else {
        #################### if there is no peak
        
        final_peakset_chr=data.frame()
        empFDR=0
    }
	
	return( list( final_peakset=final_peakset_chr, empFDR=empFDR ) )
}