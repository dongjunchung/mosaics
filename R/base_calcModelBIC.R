
.calcModelBIC <- function( loglik, n, nChr, 
	method="mosaics", analysisType="IO", signalModel="2S", type="BIC" )
{     
    # calculate number of parameters for each case
    
    npar <- switch( analysisType,
		
    	OS = 9,
    	TS = 11,
    	IO = 6
    )
    
    if ( method == "mosaicsHMM" ) {
    	npar <- npar - 1 + ( 1 + 2 ) * nChr
			# minus: proportion of binding (1)
			# plus: for each chromosome, 
			#		initial probability (1) + transition probability (2)
			# [Note] sum to one constraint
    }
    
    if ( signalModel == "2S" ) { 
    	npar <- npar + 3
    		# 1S: b, c
    		# 2S: p1, b1, c2, b2, c2
    }
    
    # calculate BIC
    
    penalty <- switch( type,
        AIC = 2,
        BIC = log(n)
    )
    
    val <- -2 * loglik + penalty * npar
    
    return(val)
}
