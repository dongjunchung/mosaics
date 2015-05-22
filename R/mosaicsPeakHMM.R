
# HMM-based MOSAiCS peak calling

setMethod(
    f="mosaicsPeakHMM",
    signature="MosaicsHMM",
    definition=function( object, FDR=0.05, decoding="posterior",
    	binsize=NA, maxgap=0, minsize=0, thres=0,
    	parallel=FALSE, nCore=8 ) 
    {     	            
        # error treatment: invalid decoding specification
        
        if ( decoding != "viterbi" && decoding != "posterior" )
        {
            stop( "Invalid decoding argument: should be either 'viterbi' or 'posterior'!" )
        }
            	            
        # parameter setting
        
        if ( decoding == "viterbi" ) {
        	message( "Info: peak calling using Viterbi algorithm." )
    	} else if ( decoding == "posterior" ) {
	    	message( "Info: peak calling using posterior decoding." )
    	}
        
        if ( is.na(binsize) ) {
            binsize <- object@binsize
        }
        
        nRatio <- object@nRatio
	    
        # HMM decoding
        
        message( "Info: calculating posterior probabilities..." )
        
        # posterior probability: P( Z_i=0 | Y )  
    
	    if ( parallel ) {
	        betapH_chr <- mclapply( object@HMMfit, 
	            function(x) {    
				    #pp0_org <- ( x$aMat[1,] * x$bMat[1,] ) / x$px
				    #pp1_org <- ( x$aMat[2,] * x$bMat[2,] ) / x$px
				    #betapH <- .ff_normalize( pp0_org, pp1_org )	
					pp0_org <- ( x$aMat[1,] * x$bMat[1,] )
					pp1_org <- ( x$aMat[2,] * x$bMat[2,] )
					betapH <- pp0_org / ( pp0_org + pp1_org )
				    return(betapH)	            
	        	}, mc.cores=nCore )    
	    } else {
	        betapH_chr <- lapply( object@HMMfit, 
	            function(x) {    
				    #pp0_org <- ( x$aMat[1,] * x$bMat[1,] ) / x$px
				    #pp1_org <- ( x$aMat[2,] * x$bMat[2,] ) / x$px
				    #betapH <- .ff_normalize( pp0_org, pp1_org )	
					pp0_org <- ( x$aMat[1,] * x$bMat[1,] )
					pp1_org <- ( x$aMat[2,] * x$bMat[2,] )
					betapH <- pp0_org / ( pp0_org + pp1_org )
				    return(betapH)	            
	        	} )
	    }
	    
	    # initial peak calling
	    
        message( "Info: calling peaks..." )
	    
	    if ( decoding == "viterbi" ) {
		    
    		# decoding using Viterbi algorithm
    		
		    if ( parallel ) {
		        bd_bin_chr <- mclapply( object@HMMfit, 
		            function(x) .ff_viterbi( x$piMat, x$gMat_chr, x$pi0Vec ),
		            mc.cores=nCore )    
		    } else {
		        bd_bin_chr <- lapply( object@HMMfit, 
		            function(x) .ff_viterbi( x$piMat, x$gMat_chr, x$pi0Vec ) 
		        	)
		    }		    
		} else if ( decoding == "posterior" ) {
			
			# decoding using posterior decoding
			
			# determine cutoff
		    	  	
		    betapH_s <- sort(unlist(betapH_chr))
		    sbetapH <- cumsum(betapH_s) / c(1:length(betapH_s))       
		    	# expected rate of false discoveries
	        id <- which( sbetapH <= FDR )
	        cutoff <- betapH_s[ max(id) ]    
		    
    		# decoding using posterior decoding
    		
	        bd_bin_chr <- lapply( betapH_chr, 
	            function(x) as.numeric( x <= cutoff )	        
	        	)		
	        
	        rm( betapH_s, sbetapH, id )
	        gc()
		}
	    
	    # polish peak lists & incorporate related information
	    
	    ann.input <- vector( "list", length(object@inputdata) )
	    
	    for ( i in 1:length(object@inputdata) ) {		    
	    	ann.input[[i]] <- list()
	    	
	    	ann.input[[i]]$inputdata <- object@inputdata[[i]]
	    	ann.input[[i]]$bd_bin <- bd_bin_chr[[i]]
	    	ann.input[[i]]$betapH <- betapH_chr[[i]]
    	}
    	names(ann.input) <- names(object@inputdata)
    
	    if ( parallel ) {
	        out <- mclapply( ann.input, 
	            function(x) .annotateHMM( 
	            	object=x, analysisType=object@peakParam@analysisType,
	            	maxgap=maxgap, minsize=minsize, thres=thres, 
	            	binsize=binsize, nRatio=nRatio ),
	            mc.cores=nCore )    
	    } else {
	        out <- lapply( ann.input, 
	            function(x) .annotateHMM( 
	            	object=x, analysisType=object@peakParam@analysisType,
	            	maxgap=maxgap, minsize=minsize, thres=thres,
	            	binsize=binsize, nRatio=nRatio )
	            )
	    }
	    
	    rm( ann.input )
	    gc()
    	
    	# post process peak lists
    	
    	peakList <- data.frame()
    	#empFDR <- 0
    	
    	for ( i in 1:length(out) ) {
	    	# peak lists
	    	
	    	if ( nrow(out[[i]]$final_peakset) > 0 ) {
		    	# stack only when this chromosome has at least one peak
		    	
		    	chr.i <- rep( names(out)[i], nrow(out[[i]]$final_peakset) )
		    	peakList.i <- data.frame( chr.i, out[[i]]$final_peakset,
		    		stringsAsFactors=FALSE )
		    	colnames(peakList.i)[1] <- "chrID"
		    	peakList <- rbind( peakList, peakList.i )
	    	
		    	rm( chr.i, peakList.i )
		    	gc()
    		}
	    }
    		
		  # binding bin, posterior probability, & empirical FDR
    	
    	chrvec <- rep( names(bd_bin_chr), sapply(bd_bin_chr,length) )
    	coordvec <- unlist(lapply( object@inputdata, function(x) x[,1] ))
    	
    	bdBin <- data.frame( chrvec, coordvec, unlist(bd_bin_chr),
    		stringsAsFactors=FALSE )
    	colnames(bdBin) <- c( "chrID", "coord", "peak" )    	
    	
    	postProb <- data.frame( chrvec, coordvec, unlist(betapH_chr),
    		stringsAsFactors=FALSE )
    	colnames(postProb) <- c( "chrID", "coord", "postProb" )
	    
    	id_bin <- which( bdBin[,3] == 1 )
	    empFDR <- sum( postProb[ id_bin, 3 ] ) / length( id_bin )
    	
    	rm( chrvec, coordvec, bd_bin_chr, betapH_chr, id_bin )
    	gc()
	        
	    message( "Info: done!" )
	    
	    # construct "MosaicsPeak" class fit	    
	            
	    peakParam <- new( "MosaicsPeakParam",
	        analysisType=object@peakParam@analysisType, 
	        signalModel=object@peakParam@signalModel, 
	        FDR=FDR, maxgap=maxgap, minsize=minsize, thres=thres,
	        decoding=decoding )
	    
      tagDataEmpty <- new( "TagData", 
          read=list(), coverage=list() )   
      
      new( "MosaicsPeak",         
	        peakList=peakList, 
          chrID=object@chrID, coord=object@coord, 
          tagCount=object@tagCount, input=object@input, 
          mappability=object@mappability, gcContent=object@gcContent,
          peakParam=peakParam, bdBin=bdBin, postProb=postProb, empFDR=empFDR,
          tagLoaded=FALSE, tagData=tagDataEmpty, seqDepth=object@seqDepth )
	}
)
