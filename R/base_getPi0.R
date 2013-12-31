
.getPi0 <- function( mu, a, Y_freq )
{
    # generate simulated tag counts
    
    b <- a/mu
    N <- length(b)    
    
    # calculate pi0
    
    if ( sum(Y_freq[ as.numeric(names(Y_freq))<=2 ])/sum(Y_freq) > 0.5 ) { 
        # if 0, 1, 2 counts are more than 50%, use 0, 1, 2 counts
        
        #print("pi0 estimation based on 0,1,2 counts")
        
        py0z0 <- dnbinom(0,a,b/(b+1)) 
        py1z0 <- dnbinom(1,a,b/(b+1))
        py2z0 <- dnbinom(2,a,b/(b+1)) 
        
        denom_pi0 <- sum(py0z0) + sum(py1z0) + sum(py2z0) 
        num_pi0 <- sum( Y_freq[1:3] )
        
        pi0 <- min( num_pi0/denom_pi0, 0.99 )
    } else {
        # otherwise, estimate pi0 by searching from mode of Y
        # until pi0 value converges or 90% bins are used
        
        # For robust estimation of modes, density estimation of ChIP tag counts is used
        
        #print("pi0 estimation based on Ymode")
    
        #Ymode <- as.numeric( names(Y_freq)[ which.max(Y_freq) ] )
        density_chip <- density( rep( as.numeric(names(Y_freq)), Y_freq ) )
        mode_smooth <- density_chip$x[ which.max(density_chip$y) ][1]
        mode_cand <- ( as.numeric(names(Y_freq))[ order( Y_freq, decreasing=TRUE ) ] )[1:10]
        Ymode <- mode_cand[ which.min( abs( mode_cand - mode_smooth ) ) ][1]
        
        if ( Ymode < 1 ) { Ymode <- as.numeric(names(Y_freq))[2] }
        
        #print( paste("mode of Y:",Ymode) )
        #print( paste("ratio Y<=Ymode:",
        #   length(which( rep( as.numeric(names(Y_freq)), Y_freq ) <= Ymode )) / sum(Y_freq) ))
        #print( paste("ratio Y>Ymode:",
        #   length(which( rep( as.numeric(names(Y_freq)), Y_freq ) > Ymode )) / sum(Y_freq) ))
        
        # estimation of pi0
        
        pi0est <- 1 # corresponding to width = -1
        cumFreq <- 0
        
        denom_pi0 <- sum( dnbinom( Ymode, a, b/(b+1) ) )
        num_pi0 <- Y_freq[ names(Y_freq)==Ymode ]
        cumFreq <- Y_freq[ names(Y_freq)==Ymode ]
                # corresponding to width = 0
         
        pi0est <- c( pi0est, num_pi0/denom_pi0 ) 
        
        #for ( i in 1:width ) {
        #   denom_pi0 <- denom_pi0 + sum( dnbinom( (Ymode-i), a, b/(b+1) ) )
        #   denom_pi0 <- denom_pi0 + sum( dnbinom( (Ymode+i), a, b/(b+1) ) )
        #   
        #   num_pi0 <- num_pi0 + Y_freq[ names(Y_freq)==(Ymode-i) ]
        #   num_pi0 <- num_pi0 + Y_freq[ names(Y_freq)==(Ymode+i) ]
        #   
        #   pi0est <- c( pi0est, num_pi0/denom_pi0 )
        #   cumFreq <- cumFreq + num_pi0
        #   
        #   print(i)
        #   print(cumFreq/sum(Y_freq))
        #} 
        
        Y_Qall <- quantile( rep( as.numeric(names(Y_freq)), Y_freq ) )
        Y_IQR <- Y_Qall[ names(Y_Qall)=="75%" ] - Y_Qall[ names(Y_Qall)=="25%" ]
        Y_Q3 <- Y_Qall[ names(Y_Qall)=="75%" ]
        
        width <- 0
        #while ( abs(pi0est[width+2]-pi0est[width+1]) > 1e-6 & cumFreq/sum(Y_freq) <= 0.90 ) {
        #while ( abs(pi0est[width+2]-pi0est[width+1]) > 1e-6 ) {
        #while ( cumFreq/sum(Y_freq) <= 0.99 ) {
        while ( width <= 2*Y_IQR ) {
            # update pi0 valuee by extending width
            # -> to avoid the case that pi0 estimate is based on too many bins
            
            width <- width + 1
            
            if ( Ymode - width > 0 ) {
                # do not use bins with tag count 0 
                
                denom_pi0 <- denom_pi0 + sum( dnbinom( (Ymode-width), a, b/(b+1) ) )
            }
            denom_pi0 <- denom_pi0 + sum( dnbinom( (Ymode+width), a, b/(b+1) ) )
            
            if ( length(which(names(Y_freq)==(Ymode-width))) > 0 & Ymode - width > 0 ) {
                # do not use bins with tag count 0 
                
                num_pi0 <- num_pi0 + Y_freq[ names(Y_freq)==(Ymode-width) ]
                cumFreq <- cumFreq + Y_freq[ names(Y_freq)==(Ymode-width) ]
            }
            if ( length(which(names(Y_freq)==(Ymode+width))) > 0 ) {
                num_pi0 <- num_pi0 + Y_freq[ names(Y_freq)==(Ymode+width) ]
                cumFreq <- cumFreq + Y_freq[ names(Y_freq)==(Ymode+width) ]
            }
            
            pi0est <- c( pi0est, num_pi0/denom_pi0 )
            
            #print(i)
            #print(num_pi0/denom_pi0)
            #print(cumFreq/sum(Y_freq))
            
            if ( abs(pi0est[width+2]-pi0est[width+1]) <= 1e-6 ) {
                # stop search if pi0 value converges
                
                #print( "pi0 value converges" )
                
                pi0 <- min( num_pi0/denom_pi0, 0.99 )
                width_final <- width
                break
            }
        }
        #if ( cumFreq/sum(Y_freq) > 0.99 ) {
        if ( width > 2*Y_IQR ) {
            # if pi0 value does not converge, then choose minimum pi0 value
            # (for robust estimation of pi0, LOESS is applied)
             
            #print( "minimum pi0 value is used" )
            
            # LOESS smoothing of pi0 estimates
            
            pi0_smooth <- loess.smooth( c((-1),0:width), pi0est )
            
            # restrict lower bound of width: Q3 - median
            # -> not trust pi0 estimate based on narrow width       
            
            pi0_smooth$y <- pi0_smooth$y[ pi0_smooth$x > (Y_Q3 - Ymode) ]
            pi0_smooth$x <- pi0_smooth$x[ pi0_smooth$x > (Y_Q3 - Ymode) ] 
            
            # determine minimum pi0 among candidate values
            
            width_final <- pi0_smooth$x[ which.min(pi0_smooth$y) ]
            pi0 <- pi0est[ which.min( abs( c((-1),0:width) - width_final ) ) ]
            pi0 <- min( pi0, 0.99 )
            #print(which.min(pi0_smooth$y))
            #print(pi0est)
        }
        #print( paste("pi0:", pi0) )
        #print( paste("width:", width_final) )
        
        #par( mfrow=c(2,1) )
        
        #hist( rep( as.numeric(names(Y_freq)), Y_freq ), nclass=10000, xlim=c(0,500) )
        #abline( v=Ymode, col="red" )
        #abline( v=Ymode-width_final, col="red", lty=2 )
        #abline( v=Ymode+width_final, col="red", lty=2 )
        #abline( v=Ymode+(Y_Q3-Ymode), col="blue" )
        #abline( v=Ymode+2*Y_IQR, col="blue" )
    
        #plot( c((-1),0:width), pi0est, type="l", xlab="Width", ylab="pi0", 
        #   main=paste("pi0:",pi0) )
        #lines( loess.smooth( c((-1),0:width), pi0est ), col="blue" )
        #abline( v=width_final, col="red" )
        #abline( v=(Y_Q3-Ymode), col="blue" )
        #abline( v=2*Y_IQR, col="blue" )
        #abline( h=pi0, col="red" )
        
        #par( mfrow=c(1,1) )
    }
    
    #pi0 <- min( num_pi0/denom_pi0, 0.99 )
    if( pi0 > 0.99 ) {
        Y_sim <- rnbinom(N,a,b/(b+1))
        pi0 <- min( num_pi0/length(which(Y_sim<=2)), 0.99 )
    }
    
    return( pi0 )
}
