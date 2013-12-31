.fitHMM <- function( inputHMM, analysisType, 
	init="mosaics", init.piMat=as.matrix(NA),
	nstate=2, binsize=200, max.iter=100, eps=1e-20 ) {
		
    # inputHMM = ( coord, peakcall, gMat0, gMat1, Y )
	
	peakcall_chr <- inputHMM[ , 2 ]
	#gMat_chr <- t(inputHMM[ , 3:4 ])
	gMat_chr <- rbind( inputHMM[ , 3 ], inputHMM[ , 4 ] )
    
    # initialization

    N <- nrow(inputHMM)
    M <- nstate
        # number of states
    aMat <- bMat <- matrix( NA, M, N )
        # forward and backward
    
    sVec <- rep( NA, N )
        # scaling factor
    
    px <- NA
    
    # initialize transition probability
        
    if ( init == "mosaics" ) {
	    # initialize based on MOSAiCS peak calling
        
	    zList <- vector( "list", N )
	    mat00 <- matrix( c( 0.7, 0.1, 0.1, 0.1 ), M, M )
	    mat01 <- matrix( c( 0.1, 0.1, 0.7, 0.1 ), M, M )
	    mat10 <- matrix( c( 0.1, 0.7, 0.1, 0.1 ), M, M )
	    mat11 <- matrix( c( 0.1, 0.1, 0.1, 0.7 ), M, M )
	    for ( i in 1:(N-1) ) {
	        if ( peakcall_chr[i] == 0 & peakcall_chr[(i+1)] == 0 ) {
	            zList[[i]] <- mat00
	        } else if ( peakcall_chr[i] == 0 & peakcall_chr[(i+1)] == 1 ) {
	            zList[[i]] <- mat01
	        } else if ( peakcall_chr[i] == 1 & peakcall_chr[(i+1)] == 0 ) {
	            zList[[i]] <- mat10
	        } else if ( peakcall_chr[i] == 1 & peakcall_chr[(i+1)] == 1 ) {
	            zList[[i]] <- mat11
	        } 
	    }
		pi0Vec <- c( 1-peakcall_chr[1], peakcall_chr[1] )  
	    piMat <- .ff_updatePiMat( zList )
    } else if ( init == "specify" ) {
	    # deterministic initialization
	    
	    pi0Vec <- c( 1-peakcall_chr[1], peakcall_chr[1] )   
	    if ( is.na(init.piMat[1,1]) ) {
	    	piMat <- matrix( c(0.9,0.1,0.1,0.9), M, M )
    	} else {
	    	piMat <- init.piMat
    	}
    } else {
	    stop( "incorrect 'init' option!" )
    }
	#print(pi0Vec)
	#print(piMat)
    
    # HMM iteration
	
	out <- .ff_emHMM( piMat, gMat_chr, pi0Vec, max.iter, eps )
    
    # calculate log likelihood
    
    loglik <- sum(log(out$sVec))
    
    # return object
	
	return( list( aMat=out$aMat, bMat=out$bMat, loglik=loglik,
		piMat=out$piMat, gMat_chr=gMat_chr, pi0Vec=out$pi0Vec ) )
}
