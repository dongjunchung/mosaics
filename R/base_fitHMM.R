.fitHMM <- function( inputHMM, peakcall, analysisType, 
	nstate=2, binsize=200, max.iter=100, eps=1e-20 ) {
		
    # extract data for given chromosome
    
    coord_chr <- inputHMM[,1]
    peakcall_chr <- inputHMM[,2]
    gMat_chr <- rbind( inputHMM[,3], inputHMM[,4] )
    Y_chr <- inputHMM[,5]
    switch( analysisType,
        OS = {
            M_chr <- inputHMM[,6]
            GC_chr <- inputHMM[,7]
        },
        TS = {
            X_chr <- inputHMM[,6]
            M_chr <- inputHMM[,7]
            GC_chr <- inputHMM[,8]
        },
        IO = {
            X_chr <- inputHMM[,6]
        }
    ) 
    
    # initialization

    N <- length(Y_chr)
    M <- nstate
        # number of states
    aMat <- bMat <- matrix( NA, M, N )
        # forward and backward
    tauMat <- matrix( 0.1, M, N )
    tauMat[ 1, peakcall_chr == 0 ] <- 0.9
    tauMat[ 2, peakcall_chr == 1 ] <- 0.9
        # membership probability
    
    pi0Vec <- c( 1-peakcall_chr[1], peakcall_chr[1] )    
        # transition from base 0
    
    piMat <- matrix( c(0.9,0.1,0.1,0.9), M, M )
        # transition
    tauList <- vector( "list", N )
    mat0 <- matrix( c(0.9,0.9,0.1,0.1), M, M )
    mat1 <- matrix( c(0.1,0.1,0.9,0.9), M, M )
    for ( i in 1:(N-1) ) {
        if ( peakcall[(i+1)] == 1 ) {
            tauList[[i]] <- mat1
        } else {
            tauList[[i]] <- mat0
        }
    }
        # transition membership probability
    
    sVec <- rep( NA, N )
        # scaling factor
    
    px <- NA
    
    # HMM iteration
    
    oldParam <- rep( 0, N )
    
    for ( iter in 1:max.iter ) {
        #print(iter)
        
        # forward
        
        fout <- .ff_forward( piMat, gMat_chr, pi0Vec )
        aMat <- fout[1:M,]
        sVec <- fout[(M+1),]        
        px <- sum( aMat[,N] )
        
        # backward
        
        bMat <- .ff_backward( piMat, gMat_chr, sVec )
        
        # E step
        
        tauList <- .ff_updateTauList( aMat, bMat, piMat, gMat_chr )        
        tauMat <- .ff_updateTauMat( tauList, gMat_chr, pi0Vec )
        
        # M step
        
        pi0Vec <- tauMat[,1]        
        piMat <- .ff_updatePiMat( tauList, tauMat )
        
        #print(pi0Vec)
        #print(piMat)
        
        # posterior probability: P( Z_j = 0 | Y )
                
        pp0_org <- ( aMat[1,] * bMat[1,] ) / px
        pp1_org <- ( aMat[2,] * bMat[2,] ) / px
        pp0 <- .ff_normalize( pp0_org, pp1_org )
        
        # stopping rule
                
        newParam <- pp0
        #print( max( abs( newParam - oldParam ) ) )
        
        if ( max( abs( newParam - oldParam ) ) < eps ) {
            break
        }
        
        oldParam <- newParam
    }
	
	return( list( aMat=aMat, bMat=bMat, px=px, 
		piMat=piMat, gMat_chr=gMat_chr, pi0Vec=pi0Vec ) )
}
