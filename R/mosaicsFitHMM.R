
# HMM-based MOSAiCS peak calling

setMethod(
    f="mosaicsFitHMM",
    signature="MosaicsFit",
    definition=function( object, signalModel="2S", binsize=NA, init="mosaics",
    	init.FDR=0.05, init.maxgap=200, init.minsize=50, init.thres=10, init.piMat=as.matrix(NA),
    	max.iter=100, eps=1e-20, parallel=FALSE, nCore=8 ) 
    {        
        init.piMat <- as.matrix(init.piMat)
		
		# check options: parallel computing (optional)
        
        if ( parallel == TRUE ) {
            message( "Use 'parallel' package for parallel computing." )        
            if ( length(find.package('parallel',quiet=TRUE)) == 0 ) {
                stop( "Please install 'parallel' package!" )
            }
        }
	        
	    # error treatment: invalid signal model specification
	    
	    if ( signalModel!="1S" && signalModel!="2S" )
	    {
	        stop( "Invalid signal model: should be either '1S' or '2S'!" )
	    }
	        
	    # error treatment: invalid initialization option
	    
	    if ( init != "mosaics" && init != "specify" )
	    {
	        stop( "Invalid initialization setting: should be either 'mosaics' or 'specify'!" )
	    }
	    
	    ##################################################################
	    #                                                                #
	    # calculation of marginal distribution: P( Y | Z )               #
	    #                                                                #
	    ##################################################################
		
		if ( init == "mosaics" ) {
			message( "Info: initialize MOSAiCS-HMM using MOSAiCS peak calling results." )
		} else if ( init == "specify" ) {
			message( "Info: initialize MOSAiCS-HMM using transition matrix user specified." )
		}
		
        
	    analysisType <- object@mosaicsEst@analysisType
	        
	    switch( signalModel,
	        "1S" = {
	            # one signal model
	            
	            message( "Info: use one-signal-component model." )
	            message( "Info: initializing HMM parameters..." )
	            
	            pi0 <- object@mosaicsEst@pi0            
	            fitMD <- .margDist_1S( mosaicsEst=object@mosaicsEst, 
	                tagCount=object@tagCount, pNfit=object@mosaicsEst@pNfit, k=3 )
	            fitPH <- .getPH_1S( margDensity=fitMD, pi0=pi0 )
	                
	            gMat <- rbind( fitMD$MDZ0, fitMD$MDZ1 )
	        },
	        "2S" = {
	            # two signal model
	            
	            message( "Info: use two-signal-component model." )
	            message( "Info: initializing HMM parameters..." )
	            
	            pi0 <- object@mosaicsEst@pi0
	            p1 <- object@mosaicsEst@p1            
	            fitMD <- .margDist_2S( mosaicsEst=object@mosaicsEst, 
	                tagCount=object@tagCount, pNfit=object@mosaicsEst@pNfit, k=3 )
	            fitPH <- .getPH_2S( margDensity=fitMD, pi0=pi0, p1=p1 )        
	                
	            gMat <- rbind( fitMD$MDZ0, fitMD$MDZ1 * p1 + fitMD$MDZ2 * (1-p1) )
	        }
	    )
	    
	    rm(fitMD)
	    gc()
	        
	    # peak calling
	    
	    #message( "Info: initializing HMM parameters..." )
	    
	    switch( analysisType,
	        OS = {
	            dataSet <- list( chrID=object@chrID, coord=object@coord, Y=object@tagCount,
	                M=object@mappability, GC=object@gcContent )
	            
	            nRatio <- 1
	        },
	        TS = {
	            dataSet <- list( chrID=object@chrID, coord=object@coord, Y=object@tagCount,
	                X=object@input, M=object@mappability, GC=object@gcContent )
	                        
              nRatio <- object@seqDepth[1] / object@seqDepth[2]
	            #nRatio <- sum(object@tagCount) / sum(object@input)
	                # genome-wide sequencing depth adjustment
	        },
	        IO = {
	            dataSet <- list( chrID=object@chrID, coord=object@coord, Y=object@tagCount, 
	                X=object@input )
	                        
	            nRatio <- object@seqDepth[1] / object@seqDepth[2]
              #nRatio <- sum(object@tagCount) / sum(object@input)
	                # genome-wide sequencing depth adjustment
	        }        
	    )
	    
	    fitPeak <- .peakCall( postProb=fitPH, dataSet=dataSet,
	        FDR=init.FDR, binsize=binsize, 
	        maxgap=init.maxgap, minsize=init.minsize, thres=init.thres, 
	        analysisType=analysisType )
	    peakcall <- fitPeak$bdBin[,3]
		#table(peakcall)
	    	# bin-level summary of peak calling results
	    
	    rm( fitPH )
	    gc()
        
	    
	    ##################################################################
	    #                                                                #
	    # HMM initialization                                             #
	    #                                                                #
	    ##################################################################
	        
	    message( "Info: estimating HMM parameters..." )       
	    
	    # split by chromosome
	    
	    #dataSet <- cbind( object@coord, peakcall, gMat[1,], gMat[2,], object@tagCount )
	    #colnames(dataSet) <- c( "coord", "peakcall", "g1", "g2", "Y" )
	    
	    switch( analysisType,
	        OS = {
	            dataSet <- cbind( object@coord, 
	            	peakcall, gMat[1,], gMat[2,], object@tagCount, 
	            	object@mappability, object@gcContent )
	            colnames(dataSet) <- c( "coord", "peakcall", "g1", "g2",
	            	"Y", "M", "GC" )
	        },
	        TS = {
	            dataSet <- cbind( object@coord,  
	            	peakcall, gMat[1,], gMat[2,], object@tagCount,
	                object@input, object@mappability, object@gcContent )
	            colnames(dataSet) <- c( "coord", "peakcall", "g1", "g2",
	            	"Y", "X", "M", "GC" )
	        },
	        IO = {
	            dataSet <- cbind( object@coord,  
	            	peakcall, gMat[1,], gMat[2,], object@tagCount, object@input )
	            colnames(dataSet) <- c( "coord", "peakcall", "g1", "g2",
	            	"Y", "X" )
	        }        
	    )
	    
	    dataSet_chr <- split( as.data.frame(dataSet), object@chrID )
	            
        # binsize calculation
        
        if ( is.na(binsize) ) {
            binsize <- min(abs( diff(object@coord) ))
        }
        
        # HMM estimation & prediction
    
	    if ( parallel ) {
	        out <- mclapply( dataSet_chr, 
	            function(x) .fitHMM( inputHMM=x, 
	            	analysisType=analysisType, init=init, init.piMat=init.piMat,
	            	nstate=2, binsize=binsize, max.iter=max.iter, eps=eps ),
	            mc.cores=nCore )    
	    } else {
	        out <- lapply( dataSet_chr, 
	            function(x) .fitHMM( inputHMM=x, 
	            	analysisType=analysisType, init=init, init.piMat=init.piMat,
	            	nstate=2, binsize=binsize, max.iter=max.iter, eps=eps )
	            )
	    }
	    
	    # BIC calculation
	    
	    message( "Info: calculating BIC of fitted models..." )
	    
    	if ( signalModel == "1S" ) {
    		BIC_mosaics <- object@bic1S
    	} else if ( signalModel == "2S" ) {
    		BIC_mosaics <- object@bic2S
    	}
    	
    	loglik <- sum( sapply( out, function(x) x$loglik ) )
    		# assume independence among chromosomes
    		# loglik = multiplication of loglik for each chromosome
    	
    	BIC_mosaicsHMM <- .calcModelBIC( 
    		loglik=loglik, n=length(object@tagCount), nChr=length(out),
    		method="mosaicsHMM", analysisType=analysisType, signalModel=signalModel, type="BIC" )
	        
	    message( "Info: done!" )
	    
	    # construct "MosaicsPeak" class fit
	            
	    peakParam <- new( "MosaicsPeakParam",
	        analysisType=analysisType, signalModel=signalModel, 
	        FDR=init.FDR, maxgap=init.maxgap, minsize=init.minsize, thres=init.thres,
	        decoding="posterior" )
	    new( "MosaicsHMM",         
	    #    HMMfit=out, mosaicsFit=object, init=init, initPiMat=init.piMat,
		      HMMfit=out, mosaicsEst=object@mosaicsEst, 
          chrID=object@chrID, coord=object@coord, 
          tagCount=object@tagCount, input=object@input, 
          mappability=object@mappability, gcContent=object@gcContent,
          inputdata=dataSet_chr, 
			    init=init, initPiMat=init.piMat,
	        peakParam=peakParam, binsize=binsize, nRatio=nRatio,
	        bicMosaics=BIC_mosaics, bicMosaicsHMM=BIC_mosaicsHMM, seqDepth=object@seqDepth )
	}
)
