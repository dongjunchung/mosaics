.decodeHMM <- function( objectHMM, decoding="viterbi", FDR=0.05, thres=0 ) {
    
    # extract fitted HMM results
    
    aMat <- objectHMM$aMat
    bMat <- objectHMM$bMat
    px <- objectHMM$px
    piMat <- objectHMM$piMat
    gMat_chr <- objectHMM$gMat_chr
    pi0Vec <- objectHMM$pi0Vec
    
    # peak calling based on HMM
	    
    # posterior probability: P( Z_i=0 | Y )  
    
    pp0_org <- ( aMat[1,] * bMat[1,] ) / px
    pp1_org <- ( aMat[2,] * bMat[2,] ) / px
    betapH <- .ff_normalize( pp0_org, pp1_org )
    
    if ( decoding == "posterior" ) {    
	    
	    # peak calling based on posterior probability
	    
	    # peak calling based on FDR
	    	  	
	    betapH_s <- sort(betapH)
	    sbetapH <- cumsum(betapH_s) / c(1:length(betapH))       # expected rate of false discoveries
        # threshold peaks by min tag count & determine peaks
        
        cutoff <- betapH_s[max(id)]    
        bd_bin <- rep( 0, length(Y) )
        bd_bin[ betapH<=cutoff & Y>thres ] <- 1
	    bd_bin <- as.numeric( sbetapH <= FDR & Y_chr >= thres )
	    
    } else if ( decoding == "viterbi" ) {
	    	    
	    # peak calling based on viterbi
	    
	    bd_bin <- .ff_viterbi( piMat, gMat_chr, pi0Vec )
	    bd_bin <- as.numeric( bd_bin == 1 & Y_chr >= thres )
	}
	
	return( list( bdBin=bd_bin, postProb=betapH ) )
}
