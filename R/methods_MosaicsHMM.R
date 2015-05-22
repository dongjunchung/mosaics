
# generic methods for "MosaicsHMM" class

setMethod(
    f="show",
    signature="MosaicsHMM",
    definition=function( object ) {    
              
        # parameters for MOSAiCS-HMM
		
		init <- object@init
		initPiMat <- object@initPiMat
		
		peakParam <- object@peakParam        
        analysisType <- peakParam@analysisType
        signalModel <- peakParam@signalModel
        FDR <- peakParam@FDR
        maxgap <- peakParam@maxgap
        minsize <- peakParam@minsize
        thres <- peakParam@thres
        
        bicMosaics <- object@bicMosaics
        bicMosaicsHMM <- object@bicMosaicsHMM
		
		# output
        
        cat( "Summary: MOSAiCS-HMM model fitting (class: MosaicsHMM)\n" )
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
            },
            "2Sr" = {
                cat( "with two robust signal components\n" )
            }
        )        
        cat( "setting for initialization:\n" )
		if ( init == "mosaics" ) {
			cat( "Initialized using MOSAiCS peak calling results.\n")
			cat( "FDR = ", FDR, ", maxgap = ", maxgap,
				", minsize = ", minsize, ", thres = ", thres, "\n", sep="" )
		} else if ( init == "specify" ) {
			cat( "initialized using the transition matrix user provided.\n")
			cat( "Transition matrix:\n" )
			cat( initPiMat[1,1], " ", initPiMat[1,2], "\n", sep="" )
			cat( initPiMat[2,1], " ", initPiMat[2,2], "\n", sep="" )
		}
        cat( "BIC of MOSAiCS model = ", bicMosaics, "\n", sep="" )
        cat( "BIC of MOSAiCS-HMM model = ", bicMosaicsHMM, "\n", sep="" )
        cat( "--------------------------------------------------\n" )
    }
)

setMethod(
    f="print",
    signature="MosaicsFit",
    definition=function( x ) {
        warning( "'print' method for 'MosaicsFit' class is not supported yet." )
    }
)

setMethod(
    f="plot",
    signature=c("MosaicsHMM","missing"),
    definition=function( x, y, seed=12345, parallel=FALSE, nCore=8, ... ) {
		# check options: parallel computing (optional)
		
		if ( parallel == TRUE ) {
			message( "Info: Use 'parallel' package for parallel computing." )        
			if ( length(find.package('parallel',quiet=TRUE)) == 0 ) {
				stop( "Please install 'parallel' package!" )
			}
		}
		
		message( "Info: summarizing ChIP-seq data..." )
	
        if ( x@mosaicsEst@analysisType=="OS" ) {
			tagCount <- unlist(lapply( x@inputdata, function(x) x$Y ))
			
            .plotGOFHMM( 
				mosaicsHMMEst=x@HMMfit, mosaicsEst=x@mosaicsEst, 
				tagCount=tagCount, input=0, 
				signalModel=x@peakParam@signalModel, k=3, 
				seed=seed, parallel=parallel, nCore=nCore )
        } else {
			tagCount <- unlist(lapply( x@inputdata, function(x) x$Y ))
			input <- unlist(lapply( x@inputdata, function(x) x$X ))
			
            .plotGOFHMM( 
				mosaicsHMMEst=x@HMMfit, mosaicsEst=x@mosaicsEst, 
				tagCount=tagCount, input=input, 
				signalModel=x@peakParam@signalModel, k=3,
				seed=seed, parallel=parallel, nCore=nCore )
        }
    }
)

setMethod(
    f="estimates",
    signature="MosaicsHMM",
    definition=function( object ) {
        mosaicsEst <- object@mosaicsEst
		
		pi0Vec <- t(sapply( object@HMMfit, function(x) x$pi0Vec ))
        piMat <- t(sapply( object@HMMfit, function(x) x$piMat ))
		piMat <- piMat[ , c(1,3,2,4) ]
		colnames(pi0Vec) <- c( "pi0", "pi1" )
		colnames(piMat) <- c( "pi00", "pi01", "pi10", "pi11" )
		
        list(            
            pi0Vec=pi0Vec,
			piMat=piMat,
            a=mosaicsEst@a,
            betaEst=mosaicsEst@betaEst,
            b=mosaicsEst@b,
            c=mosaicsEst@c,
            p1=mosaicsEst@p1,
            b1=mosaicsEst@b1,
            c1=mosaicsEst@c1,
            b2=mosaicsEst@b2,
            c2=mosaicsEst@c2
        )            
    }
)

setMethod(
    f="seqDepth",
    signature="MosaicsHMM",
    definition=function( object ) {
      return(object@seqDepth)
    }
)

